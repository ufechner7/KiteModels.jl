# apply different rel_steering values and plot turn rate
using Printf
using KiteModels, KitePodModels, KiteUtils

set = deepcopy(load_settings("system_v9.yaml"))

set.abs_tol=0.000006
set.rel_tol=0.0000001
set.elevation = 69.4
# set.v_steering = 0.2*6
# set.steering_gain = 10.0

# the following values can be changed to match your interest
set.sample_freq = 50
set.solver="DFBDF" # IDA or DFBDF
STEPS = 2400
PLOT = true
FRONT_VIEW = false
ZOOM = true
PRINT = false
STATISTIC = false
# end of user parameter section #

dt = 1/set.sample_freq

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
    last_heading = 0.0
    heading = 0.0
    for i in 1:steps
        reltime = i*dt-dt
        if reltime >= 10.0 && reltime < 10.05
            set_depower_steering(kps4.kcu, kps4.depower, steering)
        end
        last_heading = heading
        if reltime > 10.05
            az = calc_azimuth(kps4)
            heading = calc_heading(kps4)
            if heading > pi
                heading -= 2*pi
            end

            if rad2deg(heading) < -5
                set_depower_steering(kps4.kcu, kps4.depower, -steering)
            elseif rad2deg(heading) > 5
                set_depower_steering(kps4.kcu, kps4.depower, steering)
            end
        end

        KiteModels.next_step!(kps4, integrator; set_speed=0, dt)
        iter += kps4.iter
        REL_TIME[i] = reltime
        STEERING[i] = get_steering(kps4.kcu)
        HEADING[i] = rad2deg(heading)
        TURN_RATE[i] = rad2deg(heading - last_heading) / dt
        
        if plot
            if mod(i, 5) == 1
                plot2d(kps4.pos, reltime; zoom=ZOOM, front=FRONT_VIEW, segments=set.segments)                       
            end
        end
    end
end

SET_STEERING = 0.1:0.1:0.1
STEERING = zeros(STEPS)
TURN_RATE = zeros(STEPS)
HEADING = zeros(STEPS)
REL_TIME = zeros(STEPS)

for steering in 2*SET_STEERING
    integrator = KiteModels.init_sim!(kps4;  delta=0.0, stiffness_factor=1, prn=STATISTIC)
    simulate(integrator, STEPS, steering; plot=true)
end

plotx(REL_TIME, HEADING, TURN_RATE, -STEERING; xlabel="time [s]", ylabels=["heading [°]", "turn rate [deg/s]", "-rel_steering"], fig="Turn rate vs time")
# plot(REL_TIME, TURN_RATE, -STEERING; xlabel="time [s]", ylabels=["turn rate [deg/s]", "-rel_steering"], fig="Turn rate vs time")

# plot(STEERING, SIDE_CL; xlabel="rel_steering [-]", ylabel="side lift coefficient [-]", fig="Side lift coefficient vs steering")



