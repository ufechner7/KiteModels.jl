# apply different rel_steering values and plot turn rate
using Printf
using KiteModels, KitePodModels, KiteUtils

set = deepcopy(load_settings("system_v9.yaml"))

set.abs_tol=0.00006
set.rel_tol=0.000001
set.elevation = 69.4
# set.v_steering = 0.2*6
# set.steering_gain = 10.0

# the following values can be changed to match your interest
set.sample_freq = 50
set.solver="DFBDF" # IDA or DFBDF
STEPS = 3100
PLOT = true
FRONT_VIEW = true
ZOOM = true
PRINT = false
STATISTIC = false
# end of user parameter section #

dt = 1/set.sample_freq
particles = set.segments + 5
logger::Logger = Logger(particles, STEPS)

kcu::KCU = KCU(set)
kps4::KPS4 = KPS4(kcu)

if PLOT
    using Pkg
    if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
        using TestEnv; TestEnv.activate()
    end
    using ControlPlots
end

function wrap2pi(angle)
    num2pi = floor(angle / 2π + 0.5)
    angle - 2π * num2pi
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
                if rad2deg(last_heading) <= 5
                    steering +=0.1
                end
            end
        end

        KiteModels.next_step!(kps4, integrator; set_speed=0, dt)
        iter += kps4.iter
        sys_state = SysState(kps4)
        sys_state.var_16 = get_steering(kps4.kcu)
        sys_state.var_15 = rad2deg(heading - last_heading) / dt
        log!(logger, sys_state)
        
        if plot
            if mod(i, 5) == 1
                plot2d(kps4.pos, reltime; zoom=ZOOM, front=FRONT_VIEW, segments=set.segments)                       
            end
        end
    end
end

SET_STEERING = 0.1:0.1:0.1

for steering in 1*SET_STEERING
    integrator = KiteModels.init_sim!(kps4;  delta=0.0, stiffness_factor=1, prn=STATISTIC)
    simulate(integrator, STEPS, steering; plot=false)
end

function plot_steering_vs_turn_rate()
    lg = load_log("tmp")
    sl = lg.syslog
    psi = rad2deg.(wrap2pi.(sl.heading))
    plot(sl.time, -sl.var_16, sl.var_15; ylabels=["- rel_steering", "turnrate [°/s]"], fig="steering vs turnrate")
end

save_log(logger, "tmp")
if PLOT
    plot_steering_vs_turn_rate()
end
