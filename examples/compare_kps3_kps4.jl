using Printf
using KiteModels, KitePodModels, KiteUtils

set = deepcopy(load_settings("system.yaml"))
set.abs_tol=0.00006
set.rel_tol=0.000001

# the following values can be changed to match your interest
dt = 0.05
ALPHA_ZERO = 8.8            # for KPS4
STEPS = round(0.5*600/dt*0.05)
PLOT = true
FRONT_VIEW = false
ZOOM = false
PRINT = false
STATISTIC = false
# end of user parameter section #

set.version = 2

kcu::KCU  = KCU(set)
kps4::KPS4 = KPS4(kcu)
kps3::KPS3 = KPS3(kcu)

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
                plot2d(s.pos, reltime; zoom=ZOOM, xlim=(0,100), front=FRONT_VIEW, fig)          
            end
        end
    end
    iter / steps
end

integrator = KiteModels.init_sim!(kps3, delta=0.001, stiffness_factor=0.5, prn=STATISTIC)
av_steps = simulate(kps3, integrator, STEPS, true; fig="kps3")

lift, drag = KiteModels.lift_drag(kps3)
println("KPS3")
println("lift, drag  [N]: $(round(lift, digits=2)), $(round(drag, digits=2))")
println("winch_force [N]: $(round(winch_force(kps3), digits=2))")
println("Average number of callbacks per time step: $(round(av_steps, digits=2))")

kps4.set.alpha_zero = ALPHA_ZERO
integrator = KiteModels.init_sim!(kps4; delta=0.001, stiffness_factor=1, prn=STATISTIC)
av_steps = simulate(kps4, integrator, STEPS, true; fig="kps4")

lift, drag = KiteModels.lift_drag(kps4)
println("KPS4")
println("lift, drag  [N]: $(round(lift, digits=2)), $(round(drag, digits=2))")
println("winch_force [N]: $(round(winch_force(kps4), digits=2))")
println("Average number of callbacks per time step: $(round(av_steps, digits=2))")

# # KPS3
# # lift, drag  [N]: 730.25, 157.31
# # winch_force [N]: 594.89
# # Average number of callbacks per time step: 40.07
# # KPS4
# # lift, drag  [N]: 734.04, 140.91
# # winch_force [N]: 597.31
# # Average number of callbacks per time step: 101.77
