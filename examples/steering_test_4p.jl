# apply different rel_steering values and plot turn rate
using Printf
using KiteModels, KitePodModels, KiteUtils, Pkg

set = deepcopy(load_settings("system_v9.yaml"))

set.abs_tol=0.00006
set.rel_tol=0.000001
set.v_wind = 12
set.elevation = 69.4
set.l_tether = 200
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
particles = set.segments + 5
logger::Logger = Logger(particles, STEPS)

kcu::KCU = KCU(set)
kps4::KPS4 = KPS4(kcu)

if PLOT
    using Pkg
    if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
        using TestEnv; TestEnv.activate()
    end
    using ControlPlots, StatsBase
    close("all")
end

function simulate(integrator, steps; plot=false)
    OFFSET = 5
    iter = 0
    steering = 0.1
    set_depower_steering(kps4.kcu, kps4.depower, 0)
    last_heading = 0.0
    heading = 0.0
    for i in 1:steps
        reltime = i*dt-dt
        if reltime >= 10.0 && reltime < 10.05
            set_depower_steering(kps4.kcu, kps4.depower, -steering)
        end
        last_heading = heading
        if reltime > 10.05
            az = calc_azimuth(kps4)
            heading = calc_heading(kps4)
            if heading > pi
                heading -= 2*pi
            end

            if rad2deg(heading) < -OFFSET
                set_depower_steering(kps4.kcu, kps4.depower, steering)
            elseif rad2deg(heading) > OFFSET
                set_depower_steering(kps4.kcu, kps4.depower, -steering)
                if rad2deg(last_heading) <= OFFSET
                    if steering == 0.5
                        break
                    end
                    steering +=0.1
                end
            end
        end

        KiteModels.next_step!(kps4, integrator; set_speed=0, dt)
        iter += kps4.iter
        sys_state = SysState(kps4)
        sys_state.var_15 = rad2deg(heading - last_heading) / dt
        log!(logger, sys_state)
        
        if plot
            if mod(i, 5) == 1
                plot2d(kps4.pos, reltime; zoom=ZOOM, front=FRONT_VIEW, segments=set.segments)                       
            end
        end
    end
end


integrator = KiteModels.init_sim!(kps4;  delta=0.0, stiffness_factor=1, prn=STATISTIC)
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
    psi_dot = sl.var_15 # deg/s

    # p2=plot(sl.time, sl.v_app; ylabel="v_app [m/s]", fig="v_app")
    delta = delay(sl.steering, psi_dot ./ sl.v_app)
    println("delay of turnrate: $(delta*dt) s")
    delayed_steering = shift_vector(sl.steering, delta)    
    G = psi_dot ./ sl.v_app ./ delayed_steering # °/s / m/s = °/m
    for (i, g) in enumerate(G)
        if abs(delayed_steering[i]) < 0.1
            G[i] = NaN
        end
    end
    G_mean = mean(filter(!isnan, G))
    G_std = std(filter(!isnan, G))
    println("mean turnrate_law factor: $(round(G_mean, digits=3)) °/m ± $(round(G_std/G_mean*100, digits=2)) %")
    println("mean turnrate_law factor: $(round(deg2rad(G_mean), digits=4)) rad/m ± $(round(G_std/G_mean*100, digits=2)) %")
    p1 = plot(sl.time, delayed_steering, sl.var_15./sl.v_app; 
              ylabels=["delayed_steering", "turnrate/v_app [°/m]"],
              ylims=[(-0.6, 0.6), (-G_mean*0.6, G_mean*0.6)],
              fig="steering vs turnrate")
    p2 = plot(sl.time, G/G_mean; ylabel="G/G_mean [-]", fig="turnrate_law")
    display(p1); display(p2)
    return sl.time, sl.v_app, deg2rad.(psi), sl.elevation, deg2rad.(psi_dot), delayed_steering
end

function calc_c1_c2(v_app, psi, beta, psi_dot, steering)
    # c1 * v_app * u_s + c2/v_app * sin(psi) * cos(beta) - psi_dot = 0
    # Ac=b where c=[c1,c2] is a vector of the unknown coefficients.
    # A is a n×2 matrix whose columns are v_app*u_s and sin(ψ)cos(β)/vapp, and b is the vector psi_dot.
    col1 = v_app .* steering
    col2 = sin.(psi) .* cos.(beta) ./ v_app
    A = [col1 col2]
    # Then solve with c = A\b. That’ll get you a least-squares solution to the problem, 
    # the one that minimizes sum res^2, using QR decomposition.
    c = A \ psi_dot
    return c[1], c[2]
end

function plot_turnrate_law(c1, c2, time, v_app, psi, beta, psi_dot, steering)
    est_steering = psi_dot ./ (v_app * c1) .- c2 ./ (c1 .* v_app.^2) .* sin.(psi) .* cos.(beta)
    p1 = plot(time, steering, est_steering; ylabels=["delayed_steering", "est_steering"], 
              ylims=[(-0.6, 0.6), (-0.6, 0.6)],
              fig="steering vs est_steering")
    display(p1)
end

save_log(logger, "tmp")
time, v_app, psi, beta, psi_dot, steering = plot_steering_vs_turn_rate()
c1, c2 = calc_c1_c2(v_app, psi, beta, psi_dot, steering)
println("Result: c1 = $(round(c1, digits=4)) c2 = $(round(c2, digits=3))")
plot_turnrate_law(c1, c2, time, v_app, psi, beta, psi_dot, steering)

# delay of turnrate: 0.72 s
# mean turnrate_law factor: 1.4 °/m ± 20.4 %
# mean turnrate_law factor: 0.0244 rad/m ± 20.4 %
# according to Antonio: 0.048 rad/m
# according to Oriol: c1 = -0.0481175408 c2 = 1.83900976


