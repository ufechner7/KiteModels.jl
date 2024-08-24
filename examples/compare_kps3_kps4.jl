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
    if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
        using TestEnv; TestEnv.activate()
    end
    using ControlPlots
end

function simulate(s, integrator, steps, plot=false; fig="")
    iter = 0
    lines, sc, txt = nothing, nothing, nothing
    for i in 1:steps
        if PRINT
            lift, drag = KiteModels.lift_drag(s)
            @printf "%.2f: " round(integrator.t, digits=2)
            println("lift, drag  [N]: $(round(lift, digits=2)), $(round(drag, digits=2))")
        end

        KiteModels.next_step!(s, integrator; set_speed=0, dt=dt)
        iter += s.iter
        
        if plot
            reltime = i*dt-dt
            if mod(i, 5) == 1
                lines, sc, txt = plot2d(s.pos, reltime; zoom=ZOOM, front=FRONT_VIEW, segments=set.segments, fig, lines, sc, txt)    
            end
        end
    end
    iter / steps
end

integrator = KiteModels.init_sim!(kps3, delta=0, stiffness_factor=0.04, prn=STATISTIC)
av_steps = simulate(kps3, integrator, STEPS, true; fig="kps3")

lift, drag = KiteModels.lift_drag(kps3)
println("KPS3")
println("lift, drag  [N]: $(round(lift, digits=2)), $(round(drag, digits=2))")
println("winch_force [N]: $(round(winch_force(kps3), digits=2))")
println("Average number of callbacks per time step: $(round(av_steps, digits=2))")

kps4.set.alpha_zero = ALPHA_ZERO
integrator = KiteModels.init_sim!(kps4, stiffness_factor=0.5, prn=STATISTIC)
av_steps = simulate(kps4, integrator, STEPS, true; fig="kps4")

lift, drag = KiteModels.lift_drag(kps4)
println("KPS4")
println("lift, drag  [N]: $(round(lift, digits=2)), $(round(drag, digits=2))")
println("winch_force [N]: $(round(winch_force(kps4), digits=2))")
println("Average number of callbacks per time step: $(round(av_steps, digits=2))")

# KPS3
# lift, drag  [N]: 730.25, 157.31
# winch_force [N]: 594.89
# Average number of callbacks per time step: 40.07
# KPS4
# lift, drag  [N]: 733.27, 141.07
# winch_force [N]: 596.26
# Average number of callbacks per time step: 112