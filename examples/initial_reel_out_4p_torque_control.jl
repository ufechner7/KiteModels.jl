# Copyright (c) 2022, 2024 Uwe Fechner
# SPDX-License-Identifier: MIT
using Printf
using KiteModels, LinearAlgebra, SciMLBase

set = deepcopy(load_settings("system.yaml"))

# the following values can be changed to match your interest
dt = 0.010
set.solver="DFBDF" # IDA or DFBDF
set.v_reel_out = 1.0 # initial reel-out speed [m/s]
STEPS = 3000
PLOT = true
FRONT_VIEW = false
ZOOM = true
PRINT = false
STATISTIC = false
ALPHA_ZERO = 8.8 
# end of user parameter section #

set.alpha_zero = ALPHA_ZERO
set.version = 2
set.winch_model = "TorqueControlledMachine"
set.sample_freq = round(1/dt)

kcu::KCU = KCU(set)
kps4::KPS4 = KPS4(kcu)

# if PLOT
    using Pkg
    if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
        using TestEnv; TestEnv.activate()
    end
    using ControlPlots
# end

v_time = zeros(STEPS)
v_speed = zeros(STEPS)
v_set_speed = zeros(STEPS)
v_force = zeros(STEPS)
SOL = nothing

function simulate(integrator, steps, plot=false)
    global SOL
    iter = 0
    for i in 1:steps
        if PRINT
            lift, drag = KiteModels.lift_drag(kps4)
            @printf "%.2f: " round(integrator.t, digits=2)
            println("lift, drag  [N]: $(round(lift, digits=2)), $(round(drag, digits=2))")
        end
        # force needed to overcome friction
        dforce = 3.94 * set.v_reel_out
        if kps4.t_0 > 15.0
            dforce = +4.5
        end
        force = norm(kps4.forces[1])
        r = set.drum_radius
        n = set.gear_ratio
        set_torque = -r/n * force + dforce

        v_time[i] = kps4.t_0
        v_speed[i] = kps4.v_reel_out
        v_force[i] = winch_force(kps4)

        next_step!(kps4, integrator; set_torque, dt)
        if ! SciMLBase.successful_retcode(integrator.sol)
            println("Solver failed at time $(integrator.t)")
            force = norm(kps4.forces[1])
            println("force: $(force)")
            println("v_reel_out: $(kps4.v_reel_out)")
            println("sync_speed: $(kps4.sync_speed)")
            println("last_set_speed: $(kps4.wm.last_set_speed)")
            SOL = integrator.sol
            break
        end
        iter += kps4.iter
        
        if plot
            reltime = i*dt-dt
            if mod(i, 10) == 1
                plot2d(kps4.pos, reltime; zoom=ZOOM, front=FRONT_VIEW, xlim=(37, 88), 
                                        segments=set.segments, fig="side_view")            
            end
        end
    end
    iter / steps
end

kps4.sync_speed = set.v_reel_out
kps4.wm.last_set_speed = set.v_reel_out
integrator = KiteModels.init_sim!(kps4; delta=0.00015, stiffness_factor=0.8, prn=STATISTIC)

if PLOT
    av_steps = simulate(integrator, STEPS, true)
else
    println("\nStarting simulation...")
    simulate(integrator, 100)
    runtime = @elapsed av_steps = simulate(integrator, STEPS-100)
    println("\nSolver: $(set.solver)")
    println("Total simulation time: $(round(runtime, digits=3)) s")
    speed = (STEPS-100) / runtime * dt
    println("Simulation speed: $(round(speed, digits=2)) times realtime.")
end
lift, drag = KiteModels.lift_drag(kps4)
println("lift, drag  [N]: $(round(lift, digits=2)), $(round(drag, digits=2))")
println("Average number of callbacks per time step: $(round(av_steps, digits=2))")

if PLOT
    p = plotx(v_time, v_speed, v_force; 
              ylabels=["v_reelout  [m/s]","tether_force [N]"], 
              fig="winch")
    display(p)
end
# savefig("docs/src/reelout_force_4p.png")
