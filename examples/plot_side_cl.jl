# plot the side lift coefficient vs rel_steering
using Printf
using KiteModels

if haskey(ENV, "USE_V9")
    set = deepcopy(load_settings("system_v9.yaml"))
else
    set = deepcopy(load_settings("system.yaml"))
end

set.abs_tol=0.0006
set.rel_tol=0.00001
set.elevation = 69.4
set.v_steering = 0.2*6
set.steering_gain = 10.0
set.sample_freq = 50
set.depower = 38.0

# the following values can be changed to match your interest
dt = 0.02
set.solver="DFBDF" # IDA or DFBDF
STEPS = 100
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
    plt.close("all")
end

function simulate(integrator, steps, steering; plot=false)
    iter = 0
    set_depower_steering(kps4.kcu, kps4.depower, 0)
    for i in 1:steps
        if PRINT
            lift, drag = KiteModels.lift_drag(kps4)
            @printf "%.2f: " round(integrator.t, digits=2)
            println("lift, drag  [N]: $(round(lift, digits=2)), $(round(drag, digits=2))")
        end

        KiteModels.next_step!(kps4, integrator; set_speed=0, dt)
        iter += kps4.iter
        
        if plot
            reltime = i*dt-dt
            if mod(i, 5) == 1
                plot2d(kps4.pos, reltime; zoom=ZOOM, front=FRONT_VIEW, segments=set.segments)                       
            end
        end
    end
    set_depower_steering(kps4.kcu, kps4.depower, steering)
    for i in 1:15
        KiteModels.next_step!(kps4, integrator; set_speed=0, dt)
        iter += kps4.iter
    end
    kps4.side_cl, kps4.steering/kps4.set.cs_4p
end

SET_STEERING = -0.7:0.02:0.7
STEERING = zeros(length(SET_STEERING))
SIDE_CL  = zeros(length(SET_STEERING))
for (i, set_steering) in pairs(SET_STEERING)
    local side_cl, integrator
    integrator = KiteModels.init_sim!(kps4;  delta=0.001, stiffness_factor=1, prn=STATISTIC)
    side_cl, steering = simulate(integrator, STEPS, set_steering, plot=false)
    if side_cl == 0.0
        integrator = KiteModels.init_sim!(kps4;  delta=0.001, stiffness_factor=1, prn=STATISTIC)
        side_cl, steering = simulate(integrator, STEPS, set_steering, plot=false)
    end
    SIDE_CL[i] = side_cl
    STEERING[i] = steering
    println("steering: $set_steering, side_cl: $side_cl")
end
# p = plot(SET_STEERING, SIDE_CL; xlabel="rel_steering [-]", 
#          ylabel="side lift coefficient [-]", fig="Side lift coefficient vs steering")
# display(p)
p2 = plot(SET_STEERING, SIDE_CL*(set.rel_side_area/100); xlabel="set_steering [-]", 
         ylabel="side force coefficient [-]", fig="Side force coefficient vs set_steering")
display(p2)
p3 = plot(STEERING, SIDE_CL*(set.rel_side_area/100); xlabel="steering [-]", 
         ylabel="side force coefficient [-]", fig="Side force coefficient vs steering")
display(p3)


