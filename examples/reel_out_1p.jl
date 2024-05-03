using Printf
using KiteModels, KitePodModels, KiteUtils

set = deepcopy(se())

# the following values can be changed to match your interest
dt = 0.05
set.solver="DFBDF" # IDA or DFBDF
STEPS = 600
PLOT = true
FRONT_VIEW = false
ZOOM = false
PRINT = false
STATISTIC = false
# end of user parameter section #

kcu::KCU = KCU(set)
kps4::KPS4 = KPS4(kcu)
kps3::KPS3 = KPS3(kcu)

if PLOT
    using Pkg
    if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
        using TestEnv; TestEnv.activate()
    end
    using ControlPlots
    include("plot2D.jl")
end

v_time = zeros(STEPS)
v_speed = zeros(STEPS)
v_force = zeros(STEPS)

function simulate(integrator, steps, plot=false)
    start = integrator.p.iter
    line, sc, txt = nothing, nothing, nothing
    for i in 1:steps
        if PRINT
            lift, drag = KiteModels.lift_drag(kps3)
            @printf "%.2f: " round(integrator.t, digits=2)
            println("lift, drag  [N]: $(round(lift, digits=2)), $(round(drag, digits=2))")
        end
        acc = 0.0
        if kps3.t_0 > 15.0
            acc = 0.1
        end
        v_ro = kps3.sync_speed+acc*dt
        v_time[i] = kps3.t_0
        v_speed[i] = kps3.v_reel_out
        v_force[i] = winch_force(kps3)
        KiteModels.next_step!(kps3, integrator, v_ro = v_ro, dt=dt)
        if plot 
            reltime = i*dt
            if mod(i, 5) == 0
                line, sc, txt = plot2d_(kps3.pos, reltime; zoom=ZOOM, front=FRONT_VIEW, segments=se().segments, line, sc, txt)             
            end
        end
    end
    (integrator.p.iter - start) / steps
end

integrator = KiteModels.init_sim!(kps3, stiffness_factor=0.04, prn=STATISTIC)
kps3.sync_speed = 0.0

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
println("Average number of callbacks per time step: $av_steps")

# p1 = plot(v_time, v_speed, ylabel="v_reelout  [m/s]")
# p2 = plot(v_time, v_force, ylabel="tether_force [N]")
# plot(p1, p2, layout = (2, 1))
# savefig("docs/src/reelout_force_1p.png")
