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
alphas = deg2rad.(polars["alpha"])
d_flap_angles = deg2rad.(polars["d_flap_angle"])
cl_values = polars["cl"]
cd_values = polars["cd"]
c_te_values = polars["c_te"]

rm_idx = []
dist = 0.02
for i in 2:length(alphas)-1
    if d_flap_angles[i-1] == d_flap_angles[i+1] && abs(cd_values[i-1] - cd_values[i]) > dist && abs(cd_values[i+1] - cd_values[i]) > dist
        push!(rm_idx, i)
    end
end
deleteat!(alphas, rm_idx)
deleteat!(d_flap_angles, rm_idx)
deleteat!(cl_values, rm_idx)
deleteat!(cd_values, rm_idx)
deleteat!(c_te_values, rm_idx)

order = 1
println("1")
cl_spline = Spline2D(alphas, d_flap_angles, cl_values; kx = order+1, ky = order, s=1.0) # order - 2nd order s = 29 - 5th order s = 34
println("2")
cd_spline = Spline2D(alphas, d_flap_angles, cd_values; kx = order+1, ky = order, s=0.3) # order - 2nd order s = 5 - 5th order s = 8
println("3")
c_te_spline = Spline2D(alphas, d_flap_angles, c_te_values; kx = order+1, ky = order, s=0.04) # order - 2nd order s = 1 - 5th order s = 0.4

d_flap_angles_to_plot = [deg2rad(-78.5), deg2rad(2.0), deg2rad(79.5)]  # List of flap angles to plot

for d_flap_angle in d_flap_angles_to_plot
    # Filter data for the current d_flap_angle
    filtered_alphas = [alphas[i] for i in eachindex(d_flap_angles) if d_flap_angles[i] == d_flap_angle]
    @show length(filtered_alphas)
    filtered_cl_values = [cl_values[i] for i in eachindex(d_flap_angles) if d_flap_angles[i] == d_flap_angle]

    # Interpolate cl values for the filtered alphas
    spl_alphas = minimum(filtered_alphas):0.1:maximum(filtered_alphas)
    spl_cl_values = [cl_spline(alpha, d_flap_angle) for alpha in spl_alphas]

    # Plot cl_spline vs cl_values for the current d_flap_angle
    plt.plot(spl_alphas, spl_cl_values, "-", label="Interpolated Cl (Flap Angle = $d_flap_angle)")
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
    spl_alphas = minimum(filtered_alphas):0.1:maximum(filtered_alphas)
    spl_cd_values = [cd_spline(alpha, d_flap_angle) for alpha in spl_alphas]

    # Plot cd_spline vs cd_values for the current d_flap_angle
    plt.plot(spl_alphas, spl_cd_values, "-", label="Interpolated Cd (Flap Angle = $d_flap_angle)")
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
    spl_alphas = minimum(filtered_alphas):0.1:maximum(filtered_alphas)
    spl_c_te_values = [c_te_spline(alpha, d_flap_angle) for alpha in spl_alphas]

    # Plot c_te_spline vs c_te_values for the current d_flap_angle
    plt.plot(spl_alphas, spl_c_te_values, "-", label="Interpolated C_te (Flap Angle = $d_flap_angle)")
    plt.plot(filtered_alphas, filtered_c_te_values, "o", label="Actual C_te (Flap Angle = $d_flap_angle)")
end
plt.xlabel("Alpha")
plt.ylabel("C_te Values")
plt.title("Interpolated C_te vs Actual C_te for Different Flap Angles")
plt.legend()
plt.grid(true)
plt.show()

# plt.figure()
# alpha = 0.17247453081903144
# d_flap_angles = -1.0:0.01:1.0
# spl_cl_values = [cd_spline(alpha, d_flap_angle) for d_flap_angle in d_flap_angles]
# plt.plot(d_flap_angles, spl_cl_values, "-", label="Interp for different d_flap")
# plt.xlabel("Alpha")
# plt.ylabel("Cl Values")
# plt.title("Interp for different d_flap")
# plt.legend()
# plt.grid(true)
# plt.show()

@show cl_spline(-0.0866, 1.29 - 0.0866)