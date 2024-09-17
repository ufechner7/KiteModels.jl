using Distributed, Timers
using Xfoil
using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
# using ControlPlots
using KiteUtils
tic()
const procs = addprocs()

function normalize!(x, y)
    x_min = minimum(x)
    x_max = maximum(x)
    for i in eachindex(x)
        x[i] = (x[i] - x_min) / (x_max - x_min)
        y[i] = (y[i] - x_min) / (x_max - x_min)
    end
end

@everywhere begin
    using Xfoil, Statistics

    function turn_flap!(angle, x, y, lower_turn, upper_turn)
        theta = deg2rad(angle)
        x_turn = 0.75
        turn_distance = upper_turn - lower_turn
        smooth_idx = []
        rm_idx = []

        sign = theta > 0 ? 1 : -1
        y_turn = theta > 0 ? upper_turn : lower_turn
        for i in eachindex(x)
            if x_turn - turn_distance < x[i] < x_turn + turn_distance && sign * y[i] > 0
                append!(smooth_idx, i)
            elseif sign * y[i] < 0 && x_turn > x[i] > x_turn - turn_distance
                append!(rm_idx, i)
            end
            if x[i] > x_turn
                x_rel = x[i] - x_turn
                y_rel = y[i] - y_turn
                x[i] = x_turn + x_rel * cos(theta) + y_rel * sin(theta)
                y[i] = y_turn - x_rel * sin(theta) + y_rel * cos(theta)
                if theta > 0 && x[i] < x_turn - turn_distance/2 && y[i] > lower_turn
                    append!(rm_idx, i)
                elseif theta < 0 && x[i] < x_turn - turn_distance/2 && y[i] < upper_turn
                    append!(rm_idx, i)
                end
            end
        end

        #TODO: lower and upper is slightly off because of smoothing
        lower_i, upper_i = minimum(smooth_idx), maximum(smooth_idx)
        for i in smooth_idx
            window = min(i - lower_i + 1, upper_i - i + 1)
            x[i] = mean(x[i-window:i+window])
        end
        deleteat!(x, rm_idx)
        deleteat!(y, rm_idx)
        nothing
    end

    function calculate_constants(d_flap_angle_, x_, y_, cp, lower, upper)
        d_flap_angle = deg2rad(d_flap_angle_)
        x = deepcopy(x_)
        y = deepcopy(y_)
        c_te = 0.0
        if d_flap_angle > 0
            x_ref, y_ref = 0.75, upper
        else
            x_ref, y_ref = 0.75, lower
        end
        
        # straighten out the flap in order to find the trailing edge torque constant
        for i in eachindex(x)
            x_rel = x[i] - x_ref
            y_rel = y[i] - y_ref
            x[i] = x_ref + x_rel * cos(-d_flap_angle) + y_rel * sin(-d_flap_angle)
            y[i] = y_ref - x_rel * sin(-d_flap_angle) + y_rel * cos(-d_flap_angle)
        end
        
        x2 = []
        y2 = []
        for i in 1:(length(x)-1)
            if x[i] > x_ref && lower < y[i] < upper
                push!(x2, x[i])
                push!(y2, y[i])
                dx = x[i+1] - x[i]
                cp_avg = (cp[i] + cp[i+1]) / 2
                c_te -= dx * cp_avg * (x[i] - x_ref) / (1 - x_ref) # clockwise flap force at trailing edge
            end
        end
        return c_te
    end

    function solve_alpha(alphas, d_flap_angle, re, x_, y_, lower, upper, kite_speed, speed_of_sound)
        polars = []
        x = deepcopy(x_)
        y = deepcopy(y_)
        turn_flap!(d_flap_angle, x, y, lower, upper)
        Xfoil.set_coordinates(x, y)
        x, y = Xfoil.pane(npan=140)
        times_not_converged = 0
        @show d_flap_angle
        reinit = true
        for alpha in alphas
            converged = false
            cl = 0.0
            cd = 0.0
            # Solve for the given angle of attack
            cl, cd, _, _, converged = Xfoil.solve_alpha(alpha, re; iter=50, reinit=reinit, mach=kite_speed/speed_of_sound, ncrit=3)
            reinit = false
            times_not_converged += !converged
            if times_not_converged > 20
                break
            end
            if converged
                times_not_converged = 0
                _, cp = Xfoil.cpdump()
                c_te = calculate_constants(d_flap_angle, x, y, cp, lower, upper)
                push!(polars, (alpha, d_flap_angle, cl, cd, c_te))
            end
        end
        return polars
    end

    function run_solve_alpha(alphas, d_flap_angle, re, x_, y_, lower, upper, kite_speed, speed_of_sound)
        neg_alphas = sort(alphas[findall(alphas .< 0.0)], rev = true)
        pos_alphas = sort(alphas[findall(alphas .>= 0.0)])
        polars = solve_alpha(neg_alphas, d_flap_angle, re, x_, y_, lower, upper, kite_speed, speed_of_sound)
        polars = append!(polars, solve_alpha(pos_alphas, d_flap_angle, re, x_, y_, lower, upper, kite_speed, speed_of_sound))
        return polars
    end
end

function get_lower_upper(x, y)
    lower_flap = 0.0
    upper_flap = 0.0
    min_lower_distance = Inf
    min_upper_distance = Inf
    for (xi, yi) in zip(x, y)
        if yi < 0
            lower_distance = abs(xi - 0.75)
            if lower_distance < min_lower_distance
                min_lower_distance = lower_distance
                lower_flap = yi
            end
        else
            upper_distance = abs(xi - 0.75)
            if upper_distance < min_upper_distance
                min_upper_distance = upper_distance
                upper_flap = yi
            end
        end
    end
    return lower_flap, upper_flap
end

function create_polars(foil_file="naca2412.dat", polar_file="polars.csv")
    println("Creating polars")
    if !endswith(polar_file, ".csv")
        polar_file *= ".csv"
    end
    if !endswith(foil_file, ".dat")
        foil_file *= ".dat"
    end
    polar_file = joinpath(get_data_path(), polar_file)
    foil_file = joinpath(get_data_path(), foil_file)

    alphas = -180:0.5:180
    d_flap_angles = -90:0.5:90

    kite_speed = 20
    speed_of_sound = 343
    reynolds_number = kite_speed * (se("system_3l.yaml").middle_length + se("system_3l.yaml").tip_length)/2 / 1.460e-5
    println("Reynolds number for flying speed of $kite_speed is $reynolds_number")

    # Read airfoil coordinates from a file.
    x, y = open(foil_file, "r") do f
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
    Xfoil.set_coordinates(x, y)
    x, y = Xfoil.pane(npan=140)
    lower, upper = get_lower_upper(x, y)

    polars = nothing
    try
        polars = @distributed (append!) for d_flap_angle in d_flap_angles
            return run_solve_alpha(alphas, d_flap_angle, reynolds_number, x, y, lower, upper, kite_speed, speed_of_sound)
        end
    finally
        println("removing processes")
        rmprocs(procs)
    end

    println("Alpha\t\tFlap Angle\tCl\t\tCd\t\tc_te")
    for (alpha, d_flap_angle, cl, cd, c_te) in polars
        println("$alpha\t$d_flap_angle\t$(cl)\t$(cd)\t$(c_te)")
    end

    csv_content = "alpha,d_flap_angle,cl,cd,c_te\n"
    for (alpha, d_flap_angle, cl, cd, c_te) in polars
        csv_content *= string(
            alpha, ",", 
            d_flap_angle, ",", 
            cl, ",", 
            cd, ",", 
            c_te, "\n"
        )
    end
    open(polar_file, "w") do f
        write(f, csv_content)
    end
    open(foil_file, "r+") do f
        lines = readlines(f)
        if any(isletter, chomp(lines[1]))
            lines[1] *= " - polars created"
        else
            pushfirst!(lines, "polars created")
        end
        seek(f, 0)
        for line in lines
            println(f, line)
        end
    end
end

create_polars()
toc()