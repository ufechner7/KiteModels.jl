# activate the test environment if needed
using Pkg
if ! ("BenchmarkTools" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end

using Test, BenchmarkTools, StaticArrays, LinearAlgebra, KiteUtils, Plots
using KiteModels, KitePodModels

if ! @isdefined kcu
    const kcu = KCU()
end
if ! @isdefined kps
    const kps = KPS4(kcu)
end

const dt = 0.05
const PLOT = true

function plot2d(x, z; zoom=1)
    if zoom == 1
        x_max = maximum(x)
        z_max = maximum(z)
        plot(x,z, xlabel="x [m]", ylabel="z [m]", legend=false, xlims = (x_max-15.0, x_max+5), ylims = (z_max-15.0, z_max+5))
        plot!([x[7],x[10]],[z[7],z[10]], legend=false)
        plot!([x[8],x[11]],[z[8],z[11]], legend=false)
    else
        plot(x,z, xlabel="x [m]", ylabel="z [m]", legend=false)
    end
    plot!(x, z, seriestype = :scatter)
end

clear(kps)
kps.alpha_depower = deg2rad(2.2095658807330962) # from one point simulation
kps.set.alpha_zero = 0.0   
height = 134.14733504839947
kps.set.elevation = 70.7 
kps.set.profile_law = Int(EXPLOG)
kps.v_wind .= kps.v_wind_gnd * calc_wind_factor(kps, height)
kps.stiffness_factor = 1.0
# y0, yd0 = KiteModels.init_flat(kps, KiteModels.X00) # 
y0, yd0 = KiteModels.find_steady_state(kps, true)

println("\nlift, drag    [N]  : $(KiteModels.lift_drag(kps))")
println("\nSpring forces:")
spring_forces(kps)


