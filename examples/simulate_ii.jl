using Printf
using KiteModels, KitePodModels, KiteUtils

set = deepcopy(se())

set.abs_tol=0.0006
set.rel_tol=0.00001

# the following values can be changed to match your interest
dt = 0.05
STEPS = 600
PLOT = true
FRONT_VIEW = false
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
    include("plot2D.jl")
end

function simulate(integrator, steps, plot=false)
    start = integrator.p.iter
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


        KiteModels.next_step!(kps3, integrator, dt=dt)

        if plot
            reltime = i*dt
            if mod(i, 5) == 0
                lines, sc, txt = plot2d_(kps3.pos, reltime; zoom=ZOOM, front=FRONT_VIEW, segments=set.segments, lines, sc, txt)       
            end
        end
    end
    (integrator.p.iter - start) / steps
end

integrator = KiteModels.init_sim!(kps3, stiffness_factor=0.04, prn=STATISTIC)

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

# 30 to 39 times realtime in 504 iterations with integrator :Dense
# 30 to 40 times realtime in 388 iterations with integrator :GMRES