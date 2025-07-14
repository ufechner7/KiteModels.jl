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
using SymbolicAWEModels, KiteUtils, LinearAlgebra, Statistics
using ControlSystemsBase

if ! @isdefined SIMPLE
    SIMPLE = false
end

toc()

# Simulation parameters
dt = 0.05
total_time = 1.0  # Longer simulation to see oscillations
vsm_interval = 3
steps = Int(round(total_time / dt))

# Steering parameters
steering_freq = 1/2  # Hz - full left-right cycle frequency
steering_magnitude = 10.0      # Magnitude of steering input [Nm]

# Initialize model
set = Settings("system_ram.yaml")
set_values = [-50, 0.0, 0.0]  # Set values of the torques of the three winches. [Nm]

@info "Creating wing, aero, vsm_solver, sys_struct and symbolic_awe_model:"
sam = SymbolicAWEModel(set)
sam.set.abs_tol = 1e-3
sam.set.rel_tol = 1e-3
toc()

# Initialize at elevation
set.l_tethers[2] += 0.4
set.l_tethers[3] += 0.4
set.elevation = 70.0
init!(sam; remake=false, reload=false)
sys = sam.sys

@info "System initialized at:"
toc()

# Stabilize system
SymbolicAWEModels.find_steady_state!(sam; t=10.0, dt=1.0)
u0 = -sam.set.drum_radius .* sam.integrator[sys.winch_force]
sam.set_set_values(sam.integrator, u0)
simple_linearize!(sam; tstab=10.0)
lin_sam = ss(sam.simple_lin_model.A, 
             sam.simple_lin_model.B, 
             sam.simple_lin_model.C, 
             sam.simple_lin_model.D)

logger = Logger(length(sam.sys_struct.points), steps)
sys_state = SysState(sam)
t = 0.0
runtime = 0.0
integ_runtime = 0.0
bias = set.quasi_static ? 0.45 : 0.35
t0 = sam.integrator.t
set_values = zeros(3, steps)

try
    for i in 1:steps
        local steering
        global t, set_values, runtime, integ_runtime
        PLOT && plot(sam, t; zoom=false, front=false)
        
        # Calculate steering inputs based on cosine wave
        sign = t > 0.5 ? -1 : 1
        set_values[:,i] = -sam.set.drum_radius .* sam.integrator[sys.winch_force]
        set_values[:,i] .+= sign .* [10.0, steering_magnitude, -steering_magnitude]  # Opposite steering for left/right
        _vsm_interval = vsm_interval
        # Step simulation
        steptime = @elapsed next_step!(sam; set_values=set_values[:,i], 
            dt, vsm_interval=vsm_interval)
        t_new = sam.integrator.t
        integ_steptime = sam.t_step
        t = t_new - t0  # Adjust for initial stabilization time

        # Track performance after initial transient
        if (t > total_time/2)
            runtime += steptime
            integ_runtime += integ_steptime
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
save_log(logger, "tmp")
lg =load_log("tmp")
sl = lg.syslog

# Linear simulation
lin_res = lsim(lin_sam, set_values .- u0, sl.time)

# --- Updated Plotting ---
# Extract necessary data using meaningful names
force_nonlin = [sl.force[i][1] for i in eachindex(sl.force)]
len_nonlin = [[sl.l_tether[i][j] for i in eachindex(sl.l_tether)] for j in 1:3]
aoa_nonlin = rad2deg.(sl.AoA)
heading_nonlin = rad2deg.(sl.heading)
force_lin = lin_res.y[4,:] .+ force_nonlin[1]
len_lin = [lin_res.x[2+i,:] .+ len_nonlin[i][1] for i in 1:3]
aoa_lin = rad2deg.(lin_res.y[2,:]) .+ aoa_nonlin[1]
heading_lin = rad2deg.(lin_res.y[1,:]) .+ heading_nonlin[1]

p = plotx(sl.time,
    [heading_lin, heading_nonlin],
    [len_lin..., len_nonlin...],
    [aoa_lin, aoa_nonlin],
    [force_lin, force_nonlin],
    [[sl.v_app[i] - sl.v_app[1] for i in eachindex(sl.v_app)], set_values[1,:] .- set_values[1,1]];
    ylabels=["Heading [°]", "Length [m]", "AoA [°]", "Force [N]", "Params"],
    ysize=10,
    labels=[
        ["Lin", "Nonlin"],
        ["Lin1", "Lin2", "Lin3", "Nonlin1", "Nonlin2", "Nonlin3"],
        ["Lin", "Nonlin"],
        ["Lin", "Nonlin"],
        ["V_app", "Input"],
    ],
    fig="Oscillating Steering Input Response")
display(p)

@info "Performance:" times_realtime=(total_time/2)/runtime integrator_times_realtime=(total_time/2)/integ_runtime

# 55x realtime (PLOT=false, CPU: Intel i9-9980HK (16) @ 5.000GHz)
# 40-65x realtime (PLOT=false, CPU: Intel i9-9980HK (16) @ 5.000GHz) - commit 6620ed5d0a38e96930615aad9a66e4cd666955f2
# 40x realtime (PLOT=false, CPU: Intel i9-9980HK (16) @ 5.000GHz) - commit 88a78894038d3cbd50fbff83dfbe5c26266b0637
