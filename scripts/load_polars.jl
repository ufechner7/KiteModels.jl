using Dierckx
using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using ControlPlots, KiteUtils

# Load the csv file
function read_csv(filename)
    data = Dict{String, Vector{Float64}}()
    open(filename, "r") do f
        header = split(chomp(readline(f)), ",")
        for col in header
            data[col] = Float64[]
        end
        for line in eachline(f)
            values = split(chomp(line), ",")
            for (i, col) in enumerate(header)
                push!(data[col], parse(Float64, values[i]))
            end
        end
    end
    return data
end

polars = read_csv(joinpath(get_data_dir(), "polars.csv"))
alphas = polars["alpha"]
flap_angles = polars["flap_angle"]
cl_values = polars["cl"]
cd_values = polars["cd"]
c_te_values = polars["c_te"]

global smoothing_param = 1.0  # Adjust this value as needed
order = 2

global found = false
while !found
    try
        global cl_interp = Spline2D(alphas, flap_angles, cl_values; kx = order, ky = order, s=smoothing_param)
        global cd_interp = Spline2D(alphas, flap_angles, cd_values; kx = order, ky = order, s=smoothing_param)
        global c_te_interp = Spline2D(alphas, flap_angles, c_te_values; kx = order, ky = order, s=smoothing_param)
        @show smoothing_param
        global found = true
    catch e
        if isa(e, ErrorException)
            global smoothing_param += 0.1
            @show smoothing_param
        else
            rethrow(e)
        end
    end
end

flap_angles_to_plot = [-10, 0, 10, 20, 30, 45]  # List of flap angles to plot

for flap_angle in flap_angles_to_plot
    # Filter data for the current flap_angle
    filtered_alphas = [alphas[i] for i in eachindex(flap_angles) if flap_angles[i] == flap_angle]
    filtered_cl_values = [cl_values[i] for i in eachindex(flap_angles) if flap_angles[i] == flap_angle]

    # Interpolate cl values for the filtered alphas
    interpolated_alphas = minimum(filtered_alphas):0.1:maximum(filtered_alphas)
    interpolated_cl_values = [cl_interp(alpha, flap_angle) for alpha in interpolated_alphas]

    # Plot cl_interp vs cl_values for the current flap_angle
    plot(interpolated_alphas, interpolated_cl_values, "-", label="Interpolated Cl (Flap Angle = $flap_angle)")
    plot(filtered_alphas, filtered_cl_values, "o", label="Actual Cl (Flap Angle = $flap_angle)")
end

xlabel("Alpha")
ylabel("Cl Values")
title("Interpolated Cl vs Actual Cl for Different Flap Angles")
legend()
grid(true)
show()