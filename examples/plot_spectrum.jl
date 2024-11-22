

using Printf
using KiteModels, StatsBase, LinearAlgebra, DSP

if haskey(ENV, "USE_V9")
    set = deepcopy(load_settings("system_v9.yaml"))
else
    set = deepcopy(load_settings("system.yaml"))
end

using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using ControlPlots
# plt.close("all")

set.abs_tol=0.0006
set.rel_tol=0.00001
set.l_tether=200
set.v_wind = 8.0

# the following values can be changed to match your interest
dt = 0.05
fg = 2            # cut-off frequency for the filter in Hz
use_butter  = true
order = 4         # order of the Butterworth filter
set.solver="DFBDF" # IDA or DFBDF
STEPS = 600
PLOT = true
PRINT = true
STATISTIC = false
V_WIND_200    = 7.0
DEPOWER       = 0.35
F_EX = 2.23 # frequency of exertation
# end of user parameter section #

TIME = 0.0:dt:(STEPS-1)*dt

function set_tether_diameter!(se, d; c_spring_4mm = 614600, damping_4mm = 473)
    set.d_tether = d
    set.c_spring = c_spring_4mm * (d/4.0)^2
    set.damping = damping_4mm * (d/4.0)^2
end

set_tether_diameter!(set, set.d_tether)

if PLOT
    using Pkg
    if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
        using TestEnv; TestEnv.activate()
    end
    using ControlPlots
end

FORCE = zeros(STEPS)

include("filters.jl")
include("winch_controller.jl")
wcs = WinchSpeedController(dt=dt)

function simulate(kps4, integrator, logger, steps, f_ex)
    local filtered_force
    SIN = 0.5*sin.(2*π*f_ex*TIME)
    iter = 0
    last_measurement = 0.0
    butter = create_filter(fg; dt, order)
    buffer = zeros(steps)
    buffer2 = zeros(steps)
    for i in 1:steps
        force = norm(kps4.forces[1])
        FORCE[i] = force
        if use_butter
            filtered_force = apply_filter(butter, force, buffer, i)
        else
            filtered_force = ema_filter(force, last_measurement, fg, dt)
        end
        delayed_v_reelout = apply_delay(kps4.v_reel_out, buffer2, i; delay=2)
        v_set = 0.0
        set_torque = calc_set_torque(set, wcs, v_set, delayed_v_reelout, filtered_force)
        set_torque += 200*SIN[i]
        KiteModels.next_step!(kps4, integrator; set_torque, dt)
        sys_state = KiteModels.SysState(kps4)
        aoa = kps4.alpha_2
        sys_state.var_01 = aoa
        log!(logger, sys_state)
        iter += kps4.iter
    end
    nothing
end

function sim_and_plot(set; depower=DEPOWER, f_ex)
    logger = Logger(set.segments + 5, STEPS)
    set.depower = 100*depower
    set.elevation = 67.0
    kcu::KCU = KCU(set)
    kps4::KPS4 = KPS4(kcu)
    set.v_wind = V_WIND_200
    integrator = KiteModels.init_sim!(kps4; delta=0.001*0, stiffness_factor=1, prn=STATISTIC)
    set_depower_steering(kps4.kcu, depower, 0.0)
    simulate(kps4, integrator, logger, STEPS, f_ex)
    if PLOT
        p = plot(logger.time_vec, rad2deg.(logger.elevation_vec), logger.var_01_vec, xlabel="time [s]", ylabels=["elevation [°]", "aoa [°]"], 
                fig="depower: $(depower), f_ex:"*repr(f_ex))
        display(p)
        sleep(0.2)
    end
    save_log(logger, "tmp")
end

function calc_aoa_amplitude(filename)
    log = load_log(filename)
    sl  = log.syslog
    # last 4 seconds
    aoa = sl.var_01[end-(Int64(1/dt)*4):end]
    aoa = aoa .- mean(aoa)
    0.5 * (maximum(aoa) - minimum(aoa))
end

function plot_force_speed(filename)
    log = load_log(filename)
    sl  = log.syslog
    display(plot(log.syslog.time, sl.force, sl.v_reelout;
            ylabels=["force [N]", "v_reelout [m/s]"],
            fig="force_speed"*repr(set.cmq), ysize=10))
end

for f_ex in 1.5:0.1:3.0
    sim_and_plot(set; f_ex=f_ex)
    println("AOA amplitude: ", round(calc_aoa_amplitude("tmp"), digits=3), "°")
end


nothing

