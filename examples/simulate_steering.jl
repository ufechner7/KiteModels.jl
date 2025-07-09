# Copyright (c) 2022, 2024 Uwe Fechner
# SPDX-License-Identifier: MIT

using Printf
# simulate the effect of steering the kite using the one point kite model
using KiteModels

set = deepcopy(load_settings("system.yaml"))

set.abs_tol=0.0006
set.rel_tol=0.00001

# the following values can be changed to match your interest
dt = 0.05
STEPS = 600
PLOT = true
FRONT_VIEW = true
ZOOM = false
PRINT = false
STATISTIC = false
# end of user parameter section #

kcu::KCU = KCU(set)
kps3::KPS3 = KPS3(kcu)

if PLOT
    using Pkg
    if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
        using TestEnv; TestEnv.activate()
    end
    using ControlPlots
end

function simulate(integrator, steps, plot=false)
    iter = 0
    lines, sc, txt = nothing, nothing, nothing
    for i in 1:steps
        if PRINT
            lift, drag = KiteModels.lift_drag(kps3)
            @printf "%.2f: " round(integrator.t, digits=2)
            println("lift, drag  [N]: $(round(lift, digits=2)), $(round(drag, digits=2))")
        end

        if i == 300
            set_depower_steering(kps3.kcu, 0.25, 0.1)
        elseif i == 302
            set_depower_steering(kps3.kcu, 0.25, -0.1)
        elseif i == 304
            set_depower_steering(kps3.kcu, 0.25, 0.0)            
        end


        next_step!(kps3, integrator; set_speed=0, dt)
        iter += kps3.iter

        if plot
            reltime = i*dt-dt
            if mod(i, 5) == 1
                lines, sc, txt = plot2d(kps3.pos, reltime; zoom=ZOOM, front=FRONT_VIEW, segments=set.segments, 
                                        lines, sc, txt, fig="simulate_steering")       
            end
        end
    end
    iter / steps
end

integrator = KiteModels.init!(kps3; delta=0, stiffness_factor=0.04, prn=STATISTIC)

if PLOT
    av_steps = simulate(integrator, STEPS, true)
else
    println("\nStarting simulation...")
    simulate(integrator, 100)
    runtime = @elapsed av_steps = simulate(integrator, STEPS-100)
    println("\nTotal simulation time: $(round(runtime, digits=3)) s")
    speed = (STEPS-100) / runtime * dt
    println("Simulation speed: $(round(speed, digits=2)) times realtime.")
end
lift, drag = KiteModels.lift_drag(kps3)
println("lift, drag  [N]: $(round(lift, digits=2)), $(round(drag, digits=2))")
println("Average number of callbacks per time step: $(round(av_steps, digits=2))")

# 407 times realtime with 45.74 callbacks per time step