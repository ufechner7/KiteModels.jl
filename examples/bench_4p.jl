using Printf
using KiteModels, KitePodModels, KiteUtils

using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using ControlPlots

if false include("../src/KPS4.jl") end

set = deepcopy(se())

# the following values can be changed to match your interest
dt = 0.05
set.solver="DFBDF"              # IDA or DFBDF
set.linear_solver="GMRES"       # GMRES, LapackDense or Dense
STEPS = 200
PRINT = false
STATISTIC = false
PLOT=false
# end of user parameter section #

kcu::KCU = KCU(set)
kps4::KPS4 = KPS4(kcu)

v_time = zeros(STEPS)
v_speed = zeros(STEPS)
v_force = zeros(STEPS)

function simulate(integrator, steps)
    iter = 0
    for i in 1:steps
        acc = 0.0
        if kps4.t_0 > 3.0
            acc = 0.1
        end
        v_time[i] = kps4.t_0
        v_speed[i] = kps4.v_reel_out
        v_force[i] = winch_force(kps4)
        set_speed = kps4.sync_speed+acc*dt
        if PRINT
            lift, drag = KiteModels.lift_drag(kps4)
            @printf "%.2f: " round(integrator.t, digits=2)
            println("lift, drag  [N]: $(round(lift, digits=2)), $(round(drag, digits=2))")
        end

        KiteModels.next_step!(kps4, integrator; set_speed, dt=dt)
        iter += kps4.iter
    end
    iter / steps
end

integrator = KiteModels.init_sim!(kps4, stiffness_factor=0.5, prn=STATISTIC)

if PLOT
    av_steps = simulate(integrator, STEPS)
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
    p = plotx(v_time, v_speed, v_force; ylabels=["v_reelout  [m/s]","tether_force [N]"], fig="winch")
    display(p)
end

# Ryzen 7950X, Solver: DFBDF
# Total simulation time: 0.048 s
# Simulation speed: 104.06 times realtime.
# lift, drag  [N]: 454.94, 96.5
# Average number of callbacks per time step: 152.65
