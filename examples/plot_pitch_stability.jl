

using Printf
using KiteModels, LinearAlgebra, DSP

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
plt.close("all")

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
F_EX = 2.2 # frequency of exertation
# end of user parameter section #

TIME = 0.0:dt:(STEPS-1)*dt
SIN = sin.(2π*TIME*F_EX)

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

function simulate(kps4, integrator, logger, steps)
    local filtered_force
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


AOA = zeros(length(DEPOWER))
ELEV = zeros(length(DEPOWER))
V_WIND_KITE = zeros(length(DEPOWER))

elev = set.elevation
depower = DEPOWER

logger = Logger(set.segments + 5, STEPS)
set.depower = 100*depower

# set.depower_gain = 5

kcu::KCU = KCU(set)
kps4::KPS4 = KPS4(kcu)
set.elevation += 5
set.v_wind = V_WIND_200
integrator = KiteModels.init_sim!(kps4; delta=0.001*0, stiffness_factor=1, prn=STATISTIC)
if ! isnothing(integrator)
    simulate(kps4, integrator, logger, STEPS)
    set_depower_steering(kps4.kcu, depower, 0.0)
    simulate(kps4, integrator, logger, STEPS)
end
elev = rad2deg(logger.elevation_vec[end])
V_WIND_KITE= norm(kps4.v_wind)
set.elevation = 67.36

aoa = kps4.alpha_2
orient_vec = orient_euler(kps4)
alpha_depower = rad2deg(kps4.alpha_depower)
pitch = rad2deg(orient_vec[2]) + alpha_depower
v_app = norm(kps4.v_apparent)
# v_200 = calc_wind_factor(kps4.am, 200) * V_WIND
height = logger.Z_vec[end][end-2]
if PRINT
    print("Depower: $depower, alpha_dp: $(round(alpha_depower, digits=2)), aoa: $(round(aoa, digits=2)), ")
    print("pitch: $(round(pitch, digits=2))")
    println(", elevation: $(round((elev), digits=2)), height:$(round(height, digits=2))")
end
if PLOT
    p = plot(logger.time_vec, rad2deg.(logger.elevation_vec), logger.var_01_vec, xlabel="time [s]", ylabels=["elevation [°]", "aoa [°]"], 
            fig="depower: $depower")
    display(p)
    sleep(0.2)
end
save_log(logger, "tmp")

function plot_force_speed(filename)
    log = load_log(filename)
    sl  = log.syslog
    display(plot(log.syslog.time, sl.force, sl.v_reelout;
            ylabels=["force [N]", "v_reelout [m/s]"],
            fig="force_speed", ysize=10))
end

plot_force_speed("tmp")
filename = "tmp"
lg = load_log(filename)
sl = lg.syslog
nothing

