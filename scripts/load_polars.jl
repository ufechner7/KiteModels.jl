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

polars = read_csv("data/polars.csv")
alphas = polars["alpha"]
d_flap_angles = polars["d_flap_angle"]
cl_values = polars["cl"]
cd_values = polars["cd"]
c_te_values = polars["c_te"]
cl_interp = nothing
cd_interp = nothing
c_te_interp = nothing

order = 1
# factor = 10
# alphas = alphas[1:factor:end]
# d_flap_angles = d_flap_angles[1:factor:end]
# cl_values = cl_values[1:factor:end]
# cd_values = cd_values[1:factor:end]
# c_te_values = c_te_values[1:factor:end]
println("1")
cl_interp = Spline2D(alphas, d_flap_angles, cl_values; kx = order, ky = order, s=25.0) # order - 2nd order s = 29 - 5th order s = 34
println("2")
cd_interp = Spline2D(alphas, d_flap_angles, cd_values; kx = order, ky = order, s=4.0) # order - 2nd order s = 5 - 5th order s = 8
println("3")
c_te_interp = Spline2D(alphas, d_flap_angles, c_te_values; kx = order, ky = order, s=0.2) # order - 2nd order s = 1 - 5th order s = 0.4

d_flap_angles_to_plot = [-86, 2.0, 20.5]  # List of flap angles to plot

for d_flap_angle in d_flap_angles_to_plot
    # Filter data for the current d_flap_angle
    filtered_alphas = [alphas[i] for i in eachindex(d_flap_angles) if d_flap_angles[i] == d_flap_angle]
    @show length(filtered_alphas)
    filtered_cl_values = [cl_values[i] for i in eachindex(d_flap_angles) if d_flap_angles[i] == d_flap_angle]

    # Interpolate cl values for the filtered alphas
    interpolated_alphas = minimum(filtered_alphas):0.1:maximum(filtered_alphas)
    interpolated_cl_values = [cl_interp(alpha, d_flap_angle) for alpha in interpolated_alphas]

    # Plot cl_interp vs cl_values for the current d_flap_angle
    plt.plot(interpolated_alphas, interpolated_cl_values, "-", label="Interpolated Cl (Flap Angle = $d_flap_angle)")
    plt.plot(filtered_alphas, filtered_cl_values, "o", label="Actual Cl (Flap Angle = $d_flap_angle)")
end

plt.xlabel("Alpha")
plt.ylabel("Cl Values")
plt.title("Interpolated Cl vs Actual Cl for Different Flap Angles")
plt.legend()
plt.grid(true)
plt.show()

plt.figure()
for d_flap_angle in d_flap_angles_to_plot
    # Filter data for the current d_flap_angle
    filtered_alphas = [alphas[i] for i in eachindex(d_flap_angles) if d_flap_angles[i] == d_flap_angle]
    filtered_cd_values = [cd_values[i] for i in eachindex(d_flap_angles) if d_flap_angles[i] == d_flap_angle]

    # Interpolate cl values for the filtered alphas
    interpolated_alphas = minimum(filtered_alphas):0.1:maximum(filtered_alphas)
    interpolated_cd_values = [cd_interp(alpha, d_flap_angle) for alpha in interpolated_alphas]

    # Plot cd_interp vs cd_values for the current d_flap_angle
    plt.plot(interpolated_alphas, interpolated_cd_values, "-", label="Interpolated Cd (Flap Angle = $d_flap_angle)")
    plt.plot(filtered_alphas, filtered_cd_values, "o", label="Actual Cd (Flap Angle = $d_flap_angle)")
end

plt.xlabel("Alpha")
plt.ylabel("Cd Values")
plt.title("Interpolated Cd vs Actual Cd for Different Flap Angles")
plt.legend()
plt.grid(true)
plt.show()

plt.figure()
for d_flap_angle in d_flap_angles_to_plot
    # Filter data for the current d_flap_angle
    filtered_alphas = [alphas[i] for i in eachindex(d_flap_angles) if d_flap_angles[i] == d_flap_angle]
    filtered_c_te_values = [c_te_values[i] for i in eachindex(d_flap_angles) if d_flap_angles[i] == d_flap_angle]

    # Interpolate cl values for the filtered alphas
    interpolated_alphas = minimum(filtered_alphas):0.1:maximum(filtered_alphas)
    interpolated_c_te_values = [c_te_interp(alpha, d_flap_angle) for alpha in interpolated_alphas]

    # Plot c_te_interp vs c_te_values for the current d_flap_angle
    plt.plot(interpolated_alphas, interpolated_c_te_values, "-", label="Interpolated C_te (Flap Angle = $d_flap_angle)")
    plt.plot(filtered_alphas, filtered_c_te_values, "o", label="Actual C_te (Flap Angle = $d_flap_angle)")
end

plt.xlabel("Alpha")
plt.ylabel("C_te Values")
plt.title("Interpolated C_te vs Actual C_te for Different Flap Angles")
plt.legend()
plt.grid(true)
plt.show()

# @show cl_interp(rad2deg(-0.08), rad2deg(1.53))