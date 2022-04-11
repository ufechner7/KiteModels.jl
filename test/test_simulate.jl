# activate the test environment if needed
using Pkg
if ! ("BenchmarkTools" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end

using Test, BenchmarkTools, StaticArrays, LinearAlgebra, KiteUtils, Plots, Printf
using KiteModels, KitePodModels

const Model=KPS3

if ! @isdefined kcu
    const kcu = KCU()
end
if ! @isdefined kps4
    const kps4 = Model(kcu)
end

# the following values can be changed to match your interest
dt = 0.05
STEPS = 500
PLOT = true
FRONT_VIEW = true
ZOOM = true
PRINT = false
STATISTIC = false

include("plot2d.jl")

function simulate(integrator, steps, plot=false)
    start = integrator.p.iter
    for i in 1:steps
        if PRINT
            lift, drag = KiteModels.lift_drag(kps4)
            @printf "%.2f: " round(integrator.t, digits=2)
            println("lift, drag  [N]: $(round(lift, digits=2)), $(round(drag, digits=2))")
        end
        KiteModels.next_step(kps4, integrator, dt)
        if kps4.stiffness_factor < 1.0
            kps4.stiffness_factor+=0.01
        end
        if plot
            reltime = i*dt
            p = plot2d(kps4.pos, reltime; zoom=ZOOM)
            display(p)
            sleep(dt)
        end
    end
    (integrator.p.iter - start) / steps
end

integrator = KiteModels.init_sim(kps4, 1.0, STATISTIC)
kps4.stiffness_factor = 0.04

if PLOT
    simulate(integrator, 100, true)
    av_steps = simulate(integrator, STEPS-100, true)
else
    println("\nStarting simulation...")
    simulate(integrator, 100)
    runtime = @elapsed av_steps = simulate(integrator, STEPS-100)
    println("\nTotal simulation time: $(round(runtime, digits=3)) s")
    speed = (STEPS-100) / runtime * dt
    println("Simulation speed: $(round(speed, digits=2)) times realtime.")
end
println("Average number of callbacks per time step: $av_steps")
