# Copyright (c) 2024, 2025 Bart van de Lint, Uwe Fechner
# SPDX-License-Identifier: MIT

using Timers
tic()
@info "Loading packages "

PLOT = false
using Pkg
if ! ("LaTeXStrings" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using ControlPlots, LaTeXStrings
using KiteModels, LinearAlgebra, Statistics

if ! @isdefined SIMPLE
    SIMPLE = false
end

toc()

# Simulation parameters
dt = 0.05
total_time = 10.0  # Longer simulation to see oscillations
vsm_interval = 3
steps = Int(round(total_time / dt))

# Steering parameters
steering_freq = 1/2  # Hz - full left-right cycle frequency
steering_magnitude = 10.0      # Magnitude of steering input [Nm]

# Initialize model
set = Settings("system_ram.yaml")
set.segments = 3
set_values = [-50, 0.0, 0.0]  # Set values of the torques of the three winches. [Nm]
set.quasi_static = false
set.physical_model = SIMPLE ? "simple_ram" : "ram"

@info "Creating wing, aero, vsm_solver, sys_struct and symbolic_awe_model:"
sam = SymbolicAWEModel(set)
sam.set.abs_tol = 1e-2
sam.set.rel_tol = 1e-2
toc()

# Initialize at elevation
set.l_tethers[2] += 0.4
set.l_tethers[3] += 0.4
init_sim!(sam; remake=false, reload=false)
sys = sam.sys

@info "System initialized at:"
toc()

# Stabilize system
find_steady_state!(sam)

logger = Logger(length(sam.sys_struct.points), steps)
sys_state = SysState(sam)
t = 0.0
runtime = 0.0
integ_runtime = 0.0
bias = set.quasi_static ? 0.45 : 0.35
t0 = sam.integrator.t

try
    while t < total_time
        local steering
        global t, set_values, runtime, integ_runtime
        PLOT && plot(sam, t; zoom=false, front=false)
        
        # Calculate steering inputs based on cosine wave
        steering = steering_magnitude * cos(2π * steering_freq * t + bias)
        set_values = -sam.set.drum_radius .* sam.integrator[sys.winch_force]
        _vsm_interval = 1
        if t > 1.0
            set_values .+= [0.0, steering, -steering]  # Opposite steering for left/right
            _vsm_interval = vsm_interval
        end

        # Step simulation
        steptime = @elapsed next_step!(sam; set_values, dt, vsm_interval=vsm_interval)
        t_new = sam.integrator.t
        integ_steptime = sam.t_step
        t = t_new - t0  # Adjust for initial stabilization time

        # Track performance after initial transient
        if (t > total_time/2)
            runtime += steptime
            integ_runtime += integ_steptime
            sam.integrator.ps[sys.twist_damp] = 10
        end

        # Log state variables
        update_sys_state!(sys_state, sam)
        sys_state.time = t
        log!(logger, sys_state)
    end
catch e
    if isa(e, AssertionError)
        @show t
        println(e)
    else
        rethrow(e)
    end
end
@info "Total time without plotting:"
toc()

# Plot results
c = collect
save_log(logger, "tmp")
lg =load_log("tmp")
sl = lg.syslog

# --- Updated Plotting ---
# Extract necessary data using meaningful names
turn_rates_deg = rad2deg.(hcat(sl.turn_rates...))
v_reelout_23 = [sl.v_reelout[i][2] for i in eachindex(sl.v_reelout)], [sl.v_reelout[i][3] for i in eachindex(sl.v_reelout)] # Winch 2 and 3
aero_force_z = [sl.aero_force_b[i][3] for i in eachindex(sl.aero_force_b)]
aero_moment_z = [sl.aero_moment_b[i][3] for i in eachindex(sl.aero_moment_b)]
twist_angles_deg = rad2deg.(hcat(sl.twist_angles...))
AoA_deg = rad2deg.(sl.AoA)
heading_deg = rad2deg.(sl.heading)

p = plotx(sl.time,
    [turn_rates_deg[1,:], turn_rates_deg[2,:], turn_rates_deg[3,:]],
    v_reelout_23,
    [aero_force_z, aero_moment_z],
    [twist_angles_deg[1,:], twist_angles_deg[2,:], twist_angles_deg[3,:], twist_angles_deg[4,:]],
    [AoA_deg],
    [heading_deg];
    ylabels=["turn rates [°/s]", L"v_{ro}~[m/s]", "aero F/M", "twist [°]", "AoA [°]", "heading [°]"],
    ysize=10,
    labels=[
        [L"\omega_x", L"\omega_y", L"\omega_z"],
        ["v_ro[2]", "v_ro[3]"],
        [L"F_{aero,z}", L"M_{aero,z}"],
        ["twist[1]", "twist[2]", "twist[3]", "twist[4]"],
        ["AoA"],
        ["heading"]
    ],
    fig="Oscillating Steering Input Response")
display(p)

@info "Performance:" times_realtime=(total_time/2)/runtime integrator_times_realtime=(total_time/2)/integ_runtime

# 55x realtime (PLOT=false, CPU: Intel i9-9980HK (16) @ 5.000GHz)
# 40-65x realtime (PLOT=false, CPU: Intel i9-9980HK (16) @ 5.000GHz) - commit 6620ed5d0a38e96930615aad9a66e4cd666955f2
# 40x realtime (PLOT=false, CPU: Intel i9-9980HK (16) @ 5.000GHz) - commit 88a78894038d3cbd50fbff83dfbe5c26266b0637
