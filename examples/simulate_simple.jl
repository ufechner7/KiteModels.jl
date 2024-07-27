using Printf
using KiteModels, KitePodModels, KiteUtils

set = deepcopy(se())

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

function simulate(integrator, steps, plot=false)
    start = integrator.p.iter
    lines, sc, txt = nothing, nothing, nothing
    for i in 1:steps
        if PRINT
            lift, drag = KiteModels.lift_drag(kps4)
            @printf "%.2f: " round(integrator.t, digits=2)
            println("lift, drag  [N]: $(round(lift, digits=2)), $(round(drag, digits=2))")
        end

        KiteModels.next_step!(kps4, integrator; set_speed=0, dt)
        
        if plot
            reltime = i*dt-dt
            if mod(i, 5) == 1
                lines, sc, txt = plot2d(kps4.pos, reltime; zoom=ZOOM, front=FRONT_VIEW, segments=set.segments, lines, sc, txt)                       
            end
        end
    end
    (integrator.p.iter - start) / steps
end

integrator = KiteModels.init_sim!(kps4, stiffness_factor=0.5, prn=STATISTIC)

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
