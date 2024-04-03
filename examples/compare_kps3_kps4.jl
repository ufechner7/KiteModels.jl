using Printf
using KiteModels, KitePodModels, KiteUtils

set = deepcopy(se())

# the following values can be changed to match your interest
dt = 0.05
ALPHA_ZERO = 8.8 
STEPS = round(600/dt*0.05)
PLOT = true
FRONT_VIEW = false
ZOOM = true
PRINT = false
STATISTIC = false
# end of user parameter section #

kcu::KCU = KCU(set)
kps4::KPS4 = KPS4(kcu)
kps3::KPS3 = KPS3(kcu)

set.version = 2

if PLOT
    using Pkg
    if ! ("Plots" ∈ keys(Pkg.project().dependencies))
        using TestEnv; TestEnv.activate()
    end
    using Plots
    Plots.__init__()
    include("plot2d.jl")
end

function simulate(s, integrator, steps, plot=false)
    start = integrator.p.iter
    for i in 1:steps
        if PRINT
            lift, drag = KiteModels.lift_drag(s)
            @printf "%.2f: " round(integrator.t, digits=2)
            println("lift, drag  [N]: $(round(lift, digits=2)), $(round(drag, digits=2))")
        end

        KiteModels.next_step!(s, integrator, dt=dt)
        
        if plot
            reltime = i*dt
            if mod(i, 5) == 0
                p = plot2d(s.pos, reltime; zoom=ZOOM, front=FRONT_VIEW)
                display(p)                
            end
        end
    end
    (integrator.p.iter - start) / steps
end

integrator = KiteModels.init_sim!(kps3, stiffness_factor=0.04, prn=STATISTIC)
av_steps = simulate(kps3, integrator, STEPS, true)

lift, drag = KiteModels.lift_drag(kps3)
println("KPS3")
println("lift, drag  [N]: $(round(lift, digits=2)), $(round(drag, digits=2))")
println("winch_force [N]: $(round(winch_force(kps3), digits=2))")
println("Average number of callbacks per time step: $av_steps")

kps4.set.alpha_zero = ALPHA_ZERO
integrator = KiteModels.init_sim!(kps4, stiffness_factor=0.04, prn=STATISTIC)
av_steps = simulate(kps4, integrator, STEPS, true)

lift, drag = KiteModels.lift_drag(kps4)
println("KPS4")
println("lift, drag  [N]: $(round(lift, digits=2)), $(round(drag, digits=2))")
println("winch_force [N]: $(round(winch_force(kps4), digits=2))")
println("Average number of callbacks per time step: $av_steps")

# KPS3
# lift, drag  [N]: 730.25, 157.31
# winch_force [N]: 594.89
# Average number of callbacks per time step: 40.475
# KPS4
# lift, drag  [N]: 743.35, 143.21
# winch_force [N]: 607.2
# Average number of callbacks per time step: 171.185