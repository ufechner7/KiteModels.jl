# apply different rel_steering values and plot turn rate
using Printf
using KiteModels, KitePodModels, KiteUtils

set = deepcopy(load_settings("system_v9.yaml"))

set.abs_tol=0.000006
set.rel_tol=0.0000001
set.elevation = 69.4
set.v_steering = 0.2*6
# set.steering_gain = 10.0

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

function simulate(integrator, steps, steering; plot=false)
    iter = 0
    set_depower_steering(kps4.kcu, kps4.depower, 0)
    for i in 1:steps
        if PRINT
            lift, drag = KiteModels.lift_drag(kps4)
            @printf "%.2f: " round(integrator.t, digits=2)
            println("lift, drag  [N]: $(round(lift, digits=2)), $(round(drag, digits=2))")
        end
        reltime = i*dt-dt
        if reltime >= 10.0 && reltime < 10.05
            set_depower_steering(kps4.kcu, kps4.depower, steering)
        end


        KiteModels.next_step!(kps4, integrator; set_speed=0, dt)
        iter += kps4.iter
        
        if plot
            if mod(i, 5) == 1
                plot2d(kps4.pos, reltime; zoom=ZOOM, front=FRONT_VIEW, segments=set.segments)                       
            end
        end
    end
end

STEERING = 0.1:0.1:0.1
SIDE_CL = zeros(length(STEERING))

for steering in STEERING
    integrator = KiteModels.init_sim!(kps4;  delta=0.0, stiffness_factor=1, prn=STATISTIC)
    simulate(integrator, STEPS, steering; plot=true)
end

# plot(STEERING, SIDE_CL; xlabel="rel_steering [-]", ylabel="side lift coefficient [-]", fig="Side lift coefficient vs steering")



