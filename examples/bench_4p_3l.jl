using Printf
using KiteModels, KitePodModels, KiteUtils

using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using ControlPlots

set = deepcopy(load_settings("system_3l.yaml"))

# the following values can be changed to match your interest
dt = 0.05
set.solver="DFBDF" # IDA or DFBDF
STEPS = 200
PLOT = true
FRONT_VIEW = false
ZOOM = false
PRINT = false
STATISTIC = false
ALPHA_ZERO = 8.8 
# end of user parameter section #

set.alpha_zero = ALPHA_ZERO
set.version = 2

kcu::KCU = KCU(set)
kps4_3l::KPS4_3L = KPS4_3L(kcu)

v_time = zeros(STEPS)
v_speed = zeros(STEPS)
v_force = zeros(STEPS)

function simulate(integrator, steps, offset=0)
    iter = 0
    for i in 1:steps
        if PRINT
            lift, drag = KiteModels.lift_drag(kps4_3l)
            @printf "%.2f: " round(integrator.t, digits=2)
            println("lift, drag  [N]: $(round(lift, digits=2)), $(round(drag, digits=2))")
        end
        acc = 0.0
        if kps4_3l.t_0 > 3.0 + offset
            acc = 0.1
        end
        set_speeds = kps4_3l.set_speeds.+acc*dt
        v_time[i] = kps4_3l.t_0
        v_speed[i] = kps4_3l.reel_out_speeds[1]
        v_force[i] = winch_force(kps4_3l)[1]
        KiteModels.next_step!(kps4_3l, integrator; set_values=set_speeds, dt=dt)
        iter += kps4_3l.iter
    end
    println("iter: $iter", " steps: $steps")
    return iter/steps
end

integrator = KiteModels.init_sim!(kps4_3l, stiffness_factor=0.3, prn=STATISTIC)
kps4_3l.set_speeds = [0.0, 0.0, 0.0]

println("\nStarting simulation...")
simulate(integrator, 100, 100)
runtime = @elapsed av_steps = simulate(integrator, STEPS-100)
println("\nSolver: $(set.solver)")
println("Total simulation time: $(round(runtime, digits=3)) s")
speed = (STEPS-100) / runtime * dt
println("Simulation speed: $(round(speed, digits=2)) times realtime.")
if PLOT
    p = plotx(v_time[1:STEPS-100], v_speed[1:STEPS-100], v_force[1:STEPS-100]; ylabels=["v_reelout  [m/s]","tether_force [N]"], fig="winch")
    display(p)
end
lift, drag = KiteModels.lift_drag(kps4_3l)
println("lift, drag  [N]: $(round(lift, digits=2)), $(round(drag, digits=2))")
println("Average number of callbacks per time step: $av_steps")


# savefig("docs/src/reelout_force_4p.png")

# Solver: DFBDF, reltol=0.001, Ryzen desktop
# Total simulation time: 0.139 s
# Simulation speed: 36.05 times realtime.
# lift, drag  [N]: 798.46, 314.82
# Average number of callbacks per time step: 257.6
