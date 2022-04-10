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
if ! @isdefined kps4
    const kps4 = KPS4(kcu)
end

const dt = 0.05
const PLOT = true
const SHOW_FRONT = true

function plot2d(x, z; zoom=1)
    if zoom == 1
        x_max = maximum(x)
        z_max = maximum(z)
        plot(x,z, xlabel="x [m]", ylabel="z [m]", legend=false, xlims = (x_max-15.0, x_max+5), ylims = (z_max-15.0, z_max+5))
        plot!([x[7],x[10]],[z[7],z[10]], legend=false)
        plot!([x[8],x[10]],[z[8],z[10]], legend=false)
        plot!([x[8],x[11]],[z[8],z[11]], legend=false)
        plot!([x[9],x[11]],[z[9],z[11]], legend=false)
        plot!([x[9],x[7]] ,[z[9],z[7]],  legend=false)
    else
        plot(x,z, xlabel="x [m]", ylabel="z [m]", legend=false)
    end
    plot!(x, z, seriestype = :scatter)
end

function simulate(integrator, steps, plot=false)
    start = integrator.p.iter
    for i in 1:steps
        println("lift, drag    [N]  : $(KiteModels.lift_drag(kps4))")
        KiteModels.next_step(kps4, integrator, dt)
        if kps4.stiffness_factor < 1.0
            kps4.stiffness_factor+=0.01
        end
        if plot
            x = Float64[] 
            z = Float64[]
            for i in 1:length(kps4.pos)
                if SHOW_FRONT
                    push!(x, kps4.pos[i][2])
                else
                    push!(x, kps4.pos[i][1])
                end
                push!(z, kps4.pos[i][3])
            end
            p = plot2d(x, z; zoom=1)
            display(p)
            sleep(2)
        end
    end
    (integrator.p.iter - start) / steps
end

kps4.set.elevation = 70.7 
integrator = KiteModels.init_sim(kps4, 1.0)
kps4.stiffness_factor = 0.04
kps4.damping_factor = 0.5

if PLOT
    simulate(integrator, 3, true)
    av_steps = simulate(integrator, 7, true)
else
    simulate(integrator, 3)
    @time av_steps = simulate(integrator, 7)
end
println("Average number of callbacks per time step: $av_steps")
