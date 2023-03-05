using Printf
using KiteModels, KitePodModels, KiteUtils

se().abs_tol=0.00006
se().rel_tol=0.000001

if ! @isdefined kcu;  const kcu = KCU(se());   end
if ! @isdefined kps4; const kps4 = KPS4(kcu); end

# the following values can be changed to match your interest
dt = 0.05
STEPS = 600
PLOT = true
FRONT_VIEW = false
ZOOM = false
PRINT = false
STATISTIC = false
# end of user parameter section #

if PLOT
    using Pkg
    if ! ("Plots" ∈ keys(Pkg.project().dependencies))
        using TestEnv; TestEnv.activate()
    end
    using Plots
    Plots.__init__()
    include("plot2d.jl")
end

v_time = zeros(STEPS)
v_speed = zeros(STEPS)
v_force = zeros(STEPS)

function simulate(integrator, steps, plot=false)
    start = integrator.p.iter
    for i in 1:steps
        if PRINT
            lift, drag = KiteModels.lift_drag(kps4)
            @printf "%.2f: " round(integrator.t, digits=2)
            println("lift, drag  [N]: $(round(lift, digits=2)), $(round(drag, digits=2))")
        end
        acc = 0.0
        if kps4.t_0 > 15.0
            acc = 0.1
        end
        v_ro = kps4.sync_speed+acc*dt
        v_time[i] = kps4.t_0
        v_speed[i] = kps4.v_reel_out
        v_force[i] = winch_force(kps4)
        KiteModels.next_step!(kps4, integrator, v_ro = v_ro, dt=dt)
        
        if plot
            reltime = i*dt
            if mod(i, 5) == 0
                p = plot2d(kps4.pos, reltime; zoom=ZOOM, front=FRONT_VIEW, segments=se().segments)
                display(p)                
            end
        end
    end
    (integrator.p.iter - start) / steps
end

integrator = KiteModels.init_sim!(kps4, stiffness_factor=0.04, prn=STATISTIC)
kps4.sync_speed = 0.0

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
lift, drag = KiteModels.lift_drag(kps4)
println("lift, drag  [N]: $(round(lift, digits=2)), $(round(drag, digits=2))")
println("Average number of callbacks per time step: $av_steps")

p1 = plot(v_time, v_speed, ylabel="v_reelout  [m/s]", legend=false)
p2 = plot(v_time, v_force, ylabel="tether_force [N]", legend=false)
plot(p1, p2, layout = (2, 1), legend = false)
# savefig("docs/src/reelout_force_4p.png")
