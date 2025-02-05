using Printf
using KiteModels, LinearAlgebra, SciMLBase

set = deepcopy(load_settings("system.yaml"))

# the following values can be changed to match your interest
dt = 0.05
set.solver="DFBDF" # IDA or DFBDF
set.v_reel_out = 5e-10 # initial reel-out speed [m/s]
STEPS = 600
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
        dforce = 0.0
        if kps4.t_0 > 15.0
            dforce = +4.5
        end
        force = norm(kps4.forces[1])
        r = set.drum_radius
        n = set.gear_ratio
        set_torque = -r/n * force + dforce
        # println("set_torque: $set_torque")
        v_time[i] = kps4.t_0
        v_speed[i] = kps4.v_reel_out
        v_force[i] = winch_force(kps4)
        KiteModels.next_step!(kps4, integrator; set_torque, dt)
        if ! SciMLBase.successful_retcode(integrator.sol)
            println("Solver failed at time $(integrator.t)")
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
            if mod(i, 5) == 1
                plot2d(kps4.pos, reltime; zoom=ZOOM, front=FRONT_VIEW, xlim=(37, 78), 
                                        segments=set.segments, fig="side_view")            
            end
        end
    end
    iter / steps
end

# kps4.sync_speed = set.v_reel_out
# kps4.wm.last_set_speed = set.v_reel_out
integrator = KiteModels.init_sim!(kps4; delta=0, stiffness_factor=1, prn=STATISTIC)

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
    ylabels=["v_reelout  [m/s]","tether_force [N]"], fig="winch")
    display(p)
end
# savefig("docs/src/reelout_force_4p.png")

# Solver: DFBDF, reltol=0.000001
# Total simulation time: 1.049 s
# Simulation speed: 23.83 times realtime.
# lift, drag  [N]: 545.24, 102.55
# Average number of callbacks per time step: 625.604

# Solver: DFBDF, reltol=0.001
# Total simulation time: 0.178 s
# Simulation speed: 198.04 times realtime.
# lift, drag  [N]: 633.13, 119.82
# Average number of callbacks per time step: 75.67

# Solver: IDA
# Total simulation time: 0.662 s
# Simulation speed: 37.78 times realtime.
# lift, drag  [N]: 639.63, 121.71
# Average number of callbacks per time step: 379.91
