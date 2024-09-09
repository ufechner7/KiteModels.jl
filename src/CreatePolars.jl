using Distributed
using Printf
using Xfoil

const procs = addprocs()

@everywhere begin
    using Printf
    using Xfoil

    function normalize!(x, y)
        x_min = minimum(real.(x))
        x_max = maximum(real.(x))
        for i in eachindex(x)
            x[i] = (x[i] - x_min) / (x_max - x_min)
            y[i] = (y[i] - x_min) / (x_max - x_min)
        end
    end

    function turn_flap!(angle, x, y)
        theta = deg2rad(angle)
        x_turn = 0.75
        y_turn = 0.0

        # TODO: turn down around top point. then make top more smooth. opposite for downward side.

        for i in eachindex(x)
            if real(x[i]) > x_turn
                x_rel = x[i] - x_turn
                y_rel = y[i] - y_turn
                x[i] = x_turn + x_rel * cos(theta) + y_rel * sin(theta)
                y[i] = y_turn - x_rel * sin(theta) + y_rel * cos(theta)
            end
        end
        normalize!(x, y)
    end

    function calculate_constants(flap_angle_, x_, y_, cp)
        flap_angle = deg2rad(flap_angle_)
        x = deepcopy(x_)
        y = deepcopy(y_)
        c_te = 0.0
        x_ref, y_ref = 0.75, 0.0
    
        is = findall(xi -> real(xi) >= 0.75, x)
        x = x[is]
        y = y[is]
        cp = cp[is]
    
        # straighten out the flap in order to find the trailing edge torque constant
        for i in eachindex(x)
            x_rel = x[i] - x_ref
            y_rel = y[i] - y_ref
            x[i] = x_ref + x_rel * cos(-flap_angle) + y_rel * sin(-flap_angle)
            y[i] = y_ref - x_rel * sin(-flap_angle) + y_rel * cos(-flap_angle)
        end
    
        for i in 1:(length(x)-1)
            dx = x[i+1] - x[i]
            cp_avg = (cp[i] + cp[i+1]) / 2
            c_te += dx * cp_avg * (x[i] - x_ref) / (1 - x_ref) # counter-clockwise flap force at trailing edge
        end
    
        return c_te
    end

    function run_solve_alpha(alphas, d_flap_angle, re, x_, y_)
        polars = []
        x = deepcopy(x_)
        y = deepcopy(y_)
        turn_flap!(d_flap_angle, x, y)
        Xfoil.set_coordinates_cs(x, y)
        x, y = Xfoil.pane_cs(npan=140)
        times_not_converged = 0
        @show d_flap_angle
        for alpha in alphas
            converged = false
            cl = 0.0
            cd = 0.0
            # Solve for the given angle of attack
            cl, cd, _, _, converged = Xfoil.solve_alpha_cs(alpha, re; iter=50, reinit=true, mach=0.05)
            times_not_converged += !converged
            if times_not_converged > 6
                break
            end
            if converged
                times_not_converged = 0
                _, cp = Xfoil.cpdump_cs()
                c_te = calculate_constants(alpha, x, y, cp)
                flap_angle = alpha + d_flap_angle
                push!(polars, (alpha, flap_angle, cl, cd, c_te))
            end
        end
        return polars
    end
end

function create_polars(; csv_filename="polars.csv", dat_filename="naca2412.dat")
    println("Creating polars")
    if !endswith(csv_filename, ".csv")
        csv_filename *= ".csv"
    end
    if !endswith(dat_filename, ".dat")
        dat_filename *= ".dat"
    end
    if contains(csv_filename, r"[<>:\"/\\|?*]")
        error("Invalid filename: $csv_filename")
    end
    if contains(dat_filename, r"[<>:\"/\\|?*]")
        error("Invalid filename: $dat_filename")
    end
    csv_filename = joinpath(get_data_path(), csv_filename)
    dat_filename = joinpath(get_data_path(), dat_filename)

    lower = -90
    upper = 90
    step = 1
    alphas = sort(lower:step:upper, by=abs)
    d_flap_angles = sort(lower:step:upper, by=abs)
    re = 18e6  # Reynolds number

    # Read airfoil coordinates from a file.
    x, y = open(dat_filename, "r") do f
        x = Float64[]
        y = Float64[]
        for line in eachline(f)
            entries = split(chomp(line))
            try
                push!(x, parse(Float64, entries[1]))
                push!(y, parse(Float64, entries[2]))
            catch ArgumentError
            end
        end
        x, y
    end

    normalize!(x, y)

    polars = nothing
    try
        # Distribute the computation of run_solve_alpha across multiple cores
        polars = @distributed (append!) for d_flap_angle in d_flap_angles
            run_solve_alpha(alphas, d_flap_angle, re, x, y)
        end
    finally
        println("removing processes")
        rmprocs(procs)
    end

    println("Alpha\t\tFlap Angle\tCl\t\tCd\t\tc_te")
    for (alpha, flap_angle, cl, cd, c_te) in polars
        @printf("%8f\t%8f\t%8f\t%8f\t%8f\n", alpha, flap_angle, cl, cd, c_te)
    end

    csv_content = "alpha,flap_angle,cl,cd,c_te\n"
    for (alpha, flap_angle, cl, cd, c_te) in polars
        csv_content *= string(
            alpha, ",", 
            flap_angle, ",", 
            real(cl), ",", 
            real(cd), ",", 
            real(c_te), "\n"
        )
    end
    open(csv_filename, "w") do f
        write(f, csv_content)
    end
    open(dat_filename, "r+") do f
        lines = readlines(f)
        if any(isletter, chomp(lines[1]))
            println("writing")
            lines[1] *= " - polars created"
        else
            println("writing")
            pushfirst!(lines, "polars created")
        end
        seek(f, 0)
        for line in lines
            println(f, line)
        end
    end
end
