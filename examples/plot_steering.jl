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
    # println("side_force: $(kps4.side_force)")
    set_depower_steering(kps4.kcu, kps4.depower, steering)
    for i in 1:3
        KiteModels.next_step!(kps4, integrator; set_speed=0, dt)
        iter += kps4.iter
    end
    kps4.side_force[2]
end
STEERING = -0.5:0.01:0.5
SIDE_FORCE = zeros(length(STEERING))
for (i, steering) in pairs(STEERING)
    local side_force, integrator
    integrator = KiteModels.init_sim!(kps4;  delta=0.0, stiffness_factor=1, prn=STATISTIC)
    side_force = simulate(integrator, STEPS, steering, plot=false)
    SIDE_FORCE[i] = side_force
    println("steering: $steering, side_force: $side_force")
end
plot(STEERING, SIDE_FORCE; xlabel="rel_steering [-]", ylabel="side force [N]", fig="Side force vs steering")



