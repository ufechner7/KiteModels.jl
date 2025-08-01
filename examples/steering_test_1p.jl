# Copyright (c) 2022, 2024 Uwe Fechner
# SPDX-License-Identifier: MIT

# apply different rel_steering values and plot turn rate
using Printf
using KiteModels, KitePodModels, KiteUtils, Pkg

set = deepcopy(load_settings("system.yaml"))

set.abs_tol=0.00006
set.rel_tol=0.000001
set.cs_4p = 1.0
set.v_wind = 12.0
set.elevation = 69.4
set.l_tethers[1] = 200
set.depower = set.depower_offset # fully powered kite
# set.kcu_mass = 8.4
# set.v_steering = 0.2*6
# set.steering_gain = 10.0

# the following values can be changed to match your interest
set.sample_freq = 50
set.solver="DFBDF" # IDA or DFBDF
if set.kcu_model == "KCU2"
    STEPS = 2600
else
    STEPS = 2400
end
PLOT = true
FRONT_VIEW = true
ZOOM = true
PRINT = false
STATISTIC = false
# end of user parameter section #

dt = 1/set.sample_freq
particles = set.segments + 1
logger = Logger(particles, STEPS)

kcu::KCU = KCU(set)
kps3::KPS3 = KPS3(kcu)

if PLOT
    using Pkg
    if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
        using TestEnv; TestEnv.activate()
    end
    using ControlPlots, StatsBase
    close("all")
end

function simulate(integrator, steps; plot=false)
    K = 0.5
    OFFSET = 5
    iter = 0
    steering = K*0.1
    set_depower_steering(kps3.kcu, kps3.depower, 0)
    last_heading = 0.0
    heading = 0.0
    for i in 1:steps
        reltime = i*dt-dt
        if reltime >= 10.0 && reltime < 10.05
            set_depower_steering(kps3.kcu, kps3.depower, -steering)
        end
        last_heading = heading
        if reltime > 10.05
            az = calc_azimuth(kps3)
            heading = calc_heading(kps3)
            if heading > pi
                heading -= 2*pi
            end

            if rad2deg(heading) < -OFFSET
                set_depower_steering(kps3.kcu, kps3.depower, steering)
            elseif rad2deg(heading) > OFFSET
                set_depower_steering(kps3.kcu, kps3.depower, -steering)
                if rad2deg(last_heading) <= OFFSET
                    if steering ≈ K*0.5
                        break
                    end
                    steering +=K*0.1
                end
            end
        end

        next_step!(kps3, integrator; set_speed=0, dt)
        iter += kps3.iter
        sys_state = SysState(kps3)
        sys_state.var_16 = get_steering(kps3.kcu)
        sys_state.var_15 = rad2deg(heading - last_heading) / dt
        log!(logger, sys_state)
        
        if plot
            if mod(i, 5) == 1
                plot2d(kps3.pos, reltime; zoom=ZOOM, front=FRONT_VIEW, segments=set.segments,
                       fig="steering_test_1p")                       
            end
        end
    end
end

integrator = KiteModels.init!(kps3;  delta=0.0, stiffness_factor=1, prn=STATISTIC)
simulate(integrator, STEPS; plot=true)

function delay(x, y, t_max = 10)
    @assert length(x) == length(y)
    overlap = round(Int64, t_max/dt)
    z = StatsBase.crosscor(x, y, -(overlap-1):(overlap-1))
    delay_= argmax(z)
    delay_ -= overlap
    return delay_
end

# shift vector by shift to the right
function shift_vector(vec, shift)
    shift *= -1
    if shift > 0
        return [vec[1+shift:end]; zeros(shift)]
    elseif shift < 0
        return [zeros(-shift); vec[1:end+shift]]
    else
        return vec
    end
end

function plot_steering_vs_turn_rate()
    close("all")
    lg = load_log("tmp")
    sl = lg.syslog
    psi = rad2deg.(wrap2pi.(sl.heading))

    # p2=plot(sl.time, sl.v_app; ylabel="v_app [m/s]", fig="v_app")
    delta = delay(sl.var_16, sl.var_15./sl.v_app)
    println("delay of turnrate: $(delta*dt) s")
    delayed_steering = shift_vector(sl.var_16, delta)    
    G = sl.var_15./sl.v_app./delayed_steering
    for (i, g) in enumerate(G)
        if abs(delayed_steering[i]) < 0.1
            G[i] = NaN
        end
    end
    G_mean = mean(filter(!isnan, G))
    G_std = std(filter(!isnan, G))
    println("mean turnrate_law factor: $(G_mean) ± $(G_std/G_mean*100) %")
    p1 = plot(sl.time, delayed_steering, sl.var_15./sl.v_app; 
              ylabels=["delayed_steering", "turnrate/v_app [°/m]"],
              ylims=[(-0.6, 0.6), (-G_mean*0.6, G_mean*0.6)],
              fig="steering vs turnrate")
    p2 = plot(sl.time, G/G_mean; ylabel="G/G_mean [-]", fig="turnrate_law")
    display(p1); display(p2)
end

save_log(logger, "tmp")
if PLOT
    plot_steering_vs_turn_rate()
end
