using Interpolations, Statistics, Serialization, BenchmarkTools
using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using ControlPlots, KiteUtils

alphas, d_flap_angles, cl_matrix, cd_matrix, c_te_matrix = deserialize(joinpath(dirname(get_data_path()), se("system_3l.yaml").polar_file))

function replace_nan!(matrix)
    rows, cols = size(matrix)
    distance = max(rows, cols)
    for i in 1:rows
        for j in 1:cols
            if isnan(matrix[i, j])
                neighbors = []
                for d in 1:distance
                    found = false
                    if i-d >= 1 && !isnan(matrix[i-d, j]);
                        push!(neighbors, matrix[i-d, j])
                        found = true
                    end
                    if i+d <= rows && !isnan(matrix[i+d, j])
                        push!(neighbors, matrix[i+d, j])
                        found = true
                    end
                    if j-d >= 1 && !isnan(matrix[i, j-d])
                        push!(neighbors, matrix[i, j-d])
                        found = true
                    end
                    if j+d <= cols && !isnan(matrix[i, j+d])
                        push!(neighbors, matrix[i, j+d])
                        found = true
                    end
                    if found; break; end
                end
                if !isempty(neighbors)
                    matrix[i, j] = sum(neighbors) / length(neighbors)
                end
            end
        end
    end
    return nothing
end

replace_nan!(cl_matrix) # TODO: RAD2DEG
replace_nan!(cd_matrix)
replace_nan!(c_te_matrix)

cl_interp = extrapolate(scale(interpolate(cl_matrix, BSpline(Quadratic())), alphas, d_flap_angles), NaN)
cd_interp = extrapolate(scale(interpolate(cd_matrix, BSpline(Quadratic())), alphas, d_flap_angles), NaN)
c_te_interp = extrapolate(scale(interpolate(c_te_matrix, BSpline(Quadratic())), alphas, d_flap_angles), NaN)

function plot_values(alphas, d_flap_angles, matrix, interp, name)
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")

    X_data = collect(d_flap_angles) .+ zeros(length(alphas))'
    Y_data = collect(alphas)' .+ zeros(length(d_flap_angles))

    interp_matrix = similar(matrix)
    int_alphas, int_d_flap_angles = alphas .+ deg2rad(0.5), d_flap_angles .+ deg2rad(0.5)
    interp_matrix .= [interp(alpha, d_flap_angle) for alpha in int_alphas, d_flap_angle in int_d_flap_angles]
    X_int = collect(int_d_flap_angles) .+ zeros(length(int_alphas))'
    Y_int = collect(int_alphas)' .+ zeros(length(int_d_flap_angles))

    ax.plot_wireframe(X_data, Y_data, matrix, edgecolor="royalblue", lw=0.5, rstride=5, cstride=5, alpha=0.6)
    ax.plot_wireframe(X_int, Y_int, interp_matrix, edgecolor="orange", lw=0.5, rstride=5, cstride=5, alpha=0.6)
    plt.xlabel("Alpha")
    plt.ylabel("Flap angle")
    plt.zlabel("$name values")
    plt.title("$name for different d_flap and angle")
    plt.legend()
    plt.grid(true)
    plt.show()
end


plot_values(alphas, d_flap_angles, cl_matrix, cl_interp, "Cl")
plot_values(alphas, d_flap_angles, cd_matrix, cd_interp, "Cd")
plot_values(alphas, d_flap_angles, c_te_matrix, c_te_interp, "C_te")
display(plot(alphas, cd_interp.(alphas, 0.0)))
# @show gradient(cd_interp, 10, 0)

# @benchmark cd_interp(rad2deg(rand()),rad2deg(rand()))
# Dierckx
# @benchmark cd_interp(rand(),rand())
# Range (min … max):  156.364 ns … 102.205 μs  ┊ GC (min … max): 0.00% … 99.77%
# Time  (median):     218.206 ns               ┊ GC (median):    0.00%
# Time  (mean ± σ):   223.494 ns ±   1.021 μs  ┊ GC (mean ± σ):  5.06% ±  2.26%

#                ▅                   █▅ ▃▃                        
#  ▂▁▂▂▂▂▂▂▂▂▂▃▃▂█▄▄▅▂▂▂▂▂▂▂▂▂▂▂▂▂▃▃▃██▇██▆▄▃▂▂▂▂▂▂▂▂▂▂▂▂▂▂▁▂▁▂▂ ▃
#  156 ns           Histogram: frequency by time          264 ns <
# Memory estimate: 160 bytes, allocs estimate: 4.

# Interpolations.jl
# BenchmarkTools.Trial: 10000 samples with 985 evaluations.
#  Range (min … max):  53.683 ns … 271.883 ns  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     69.062 ns               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   73.499 ns ±  15.021 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%

#            ▄▇██▇▅▅▃▃▂▁▁▁▁▁                             ▁▂▁▁    ▂
#   ▇▇█▇▅▄▃▂███████████████████▇▆▇▇▆▅▆▅▆▆▆▆▆▅▄▆▆▇▆▆▆▇▆▆▆█████▇█▆ █
#   53.7 ns       Histogram: log(frequency) by time       128 ns <
#  Memory estimate: 48 bytes, allocs estimate: 3.