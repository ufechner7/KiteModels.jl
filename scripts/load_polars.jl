using Dierckx, Statistics
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

wd = 2.0 .- abs.(cd_values ./ argmax(abs, cd_values))
order = 2
println("1")
cl_spline = Spline2D(alphas, d_flap_angles, cl_values; kx=order, ky=order, s=20.0)
println("2")
cd_spline = Spline2D(alphas, d_flap_angles, cd_values; w=wd, kx=order, ky=order, s=10.0)
println("3")
c_te_spline = Spline2D(alphas, d_flap_angles, c_te_values; kx=order, ky=order, s=1.0)

function plot_values(spline, values, name)
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    # alpha = 0.6396476065600847
    plot_alphas = 0.0:0.04:deg2rad(60)
    plot_flap_angles = -1.5:0.1:deg2rad(60)

    idx = reduce(vcat, [findall(x -> abs(x - alpha) < deg2rad(0.3), alphas) for alpha in plot_alphas])
    
    spl_values = reduce(vcat, [reduce(vcat, [spline(alpha, flap_angle) for flap_angle in plot_flap_angles]) for alpha in plot_alphas])
    extended_plot_alphas = reduce(vcat, [reduce(vcat, [alpha for _ in plot_flap_angles]) for alpha in plot_alphas])
    extended_plot_flap_angles = reduce(vcat, [reduce(vcat, [flap_angle for flap_angle in plot_flap_angles]) for _ in plot_alphas])

    ax.scatter(d_flap_angles[idx], alphas[idx], values[idx])
    ax.scatter(extended_plot_flap_angles, extended_plot_alphas, spl_values)
    plt.xlabel("Flap angle")
    plt.ylabel("Alpha")
    plt.zlabel("$name values")
    plt.title("$name for different d_flap and angle")
    plt.legend()
    plt.grid(true)
    plt.show()
end

plot_values(cd_spline, cd_values, "Cd")
plot_values(cl_spline, cl_values, "Cl")
plot_values(c_te_spline, c_te_values, "C_te")

cd_spline(0,0)*π