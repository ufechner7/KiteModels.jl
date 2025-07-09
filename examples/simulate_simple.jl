# Copyright (c) 2022, 2024 Uwe Fechner
# SPDX-License-Identifier: MIT

using Printf
# simple parking test without changing the control input
# shows how to log, plot, and print the simulation results

using KiteModels, LinearAlgebra

set = deepcopy(load_settings("system.yaml"))

set.abs_tol=0.0006
set.rel_tol=0.00001

# the following values can be changed to match your interest
dt = 0.05
set.solver="DFBDF" # IDA or DFBDF
STEPS = 600
PLOT = true
FRONT_VIEW = false
ZOOM = true
PRINT = false
STATISTIC = false
UPWIND_DIR = -pi/2 +deg2rad(10)
# end of user parameter section #

kcu::KCU = KCU(set)
kps4::KPS4 = KPS4(kcu)

if PLOT
    using Pkg
    if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
        using TestEnv; TestEnv.activate()
    end
    using ControlPlots
end

logger = Logger(set.segments + 5, STEPS)

function simulate(integrator, steps, plot=false)
    iter = 0
    for i in 1:steps
        if PRINT
            lift, drag = KiteModels.lift_drag(kps4)
            @printf "%.2f: " round(integrator.t, digits=2)
            println("lift, drag  [N]: $(round(lift, digits=2)), $(round(drag, digits=2))")
        end

        next_step!(kps4, integrator; set_speed=0, upwind_dir=UPWIND_DIR, dt)
        iter += kps4.iter    
        if plot
            reltime = i*dt-dt
            if mod(i, 5) == 1
                plot2d(kps4.pos, reltime; zoom=ZOOM, xlim=(40,60), front=FRONT_VIEW, 
                segments=set.segments, fig="upwind_dir = $(rad2deg(UPWIND_DIR)) °")                       
            end
        end
        sys_state = SysState(kps4)
        sys_state.var_01 = kps4.pitch
        sys_state.var_02 = kps4.pitch_rate
        log!(logger, sys_state)
    end
    iter / steps
end
kps4.set.upwind_dir = rad2deg(UPWIND_DIR)
integrator = KiteModels.init_sim!(kps4;  delta=0.0, stiffness_factor=1, prn=STATISTIC)

if PLOT
    global flight_log
    av_steps = simulate(integrator, STEPS, true)
    flight_log = KiteUtils.sys_log(logger)
    p = plotx(flight_log.syslog.time, flight_log.z, 
              rad2deg.(flight_log.syslog.var_01), rad2deg.(flight_log.syslog.var_02);
              xlabel="time [s]", ylabels=["z [m]", "pitch [°]", "pitch_rate [°/s]"], 
              fig="plot_height_pitch")
    display(p)
else
    println("\nStarting simulation...")
    simulate(integrator, 100)
    runtime = @elapsed av_steps = simulate(integrator, STEPS-100)
    println("\nTotal simulation time: $(round(runtime, digits=3)) s")
    speed = (STEPS-100) / runtime * dt
    println("Simulation speed: $(round(speed, digits=2)) times realtime.")
end
lift, drag = KiteModels.lift_drag(kps4)
println("Ground wind speed: $(round(set.v_wind, digits=2)) m/s")
println("lift, drag  [N]: $(round(lift, digits=2)), $(round(drag, digits=2))")
println("Average number of callbacks per time step: $(round(av_steps, digits=2))")

points = KiteUtils.get_particles(set.height_k, set.h_bridle, set.width, set.m_k, [0, 0, 0], [0, 0, -1], [10, 0, 0])
pos_A = points[3]
pos_C = points[5]
pos_D = points[6]
Pc = 0.5*(pos_C + pos_D)
distance = norm(pos_A-Pc)
println("Distance between A and Pc: $(round(distance, digits=2)) m")
