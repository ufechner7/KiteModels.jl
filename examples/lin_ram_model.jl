# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

#=
This example demonstrates linearized model accuracy by comparing:
1. Nonlinear RamAirKite model simulation 
2. Linearized state-space model simulation

Both models start from the same operating point and are subjected
to identical steering inputs. The resulting state trajectories are
plotted together to visualize how well the linearized model
approximates the nonlinear dynamics.
=#

using Timers
tic()
@info "Loading packages "

using KiteModels, LinearAlgebra, Statistics, OrdinaryDiffEqCore 
using ModelingToolkit

PLOT = true
if PLOT
    using Pkg
    if ! ("LaTeXStrings" ∈ keys(Pkg.project().dependencies))
        using TestEnv; TestEnv.activate()
    end
    using ControlPlots, LaTeXStrings
    import ControlPlots: plot
end
toc()

include(joinpath(@__DIR__, "plotting.jl"))

# TODO: use sparse autodiff

# Simulation parameters
dt = 0.001
total_time = 1.0  # Increased from 0.1s to 1.0s for better dynamics observation
vsm_interval = 3
steps = Int(round(total_time / dt))

# Steering parameters
steering_freq = 1/2  # Hz - full left-right cycle frequency
steering_magnitude = 5.0      # Magnitude of steering input [Nm]

# Initialize model
set = load_settings("system_ram.yaml")
set.segments = 3
set_values = [-50.0, 0.0, 0.0]  # Set values of the torques of the three winches. [Nm]
set.quasi_static = false
set.physical_model = "simple_ram"

@info "Creating RamAirKite model..."
s = RamAirKite(set)
s.set.abs_tol = 1e-2
s.set.rel_tol = 1e-2
toc()

measure = Measurement()

# Define outputs for linearization - angular velocities
@variables ω_b(ModelingToolkit.t_nounits)[1:3]

# Initialize at elevation with linearization outputs
s.point_system.winches[2].tether_length += 0.2
s.point_system.winches[3].tether_length += 0.2
measure.sphere_pos .= deg2rad.([65.0 65.0; 1.0 -1.0])
KiteModels.init_sim!(s, measure; 
    remake=false,
    reload=true,
    lin_outputs=[ω_b...]  # Specify which outputs to track in linear model
)
sys = s.sys

@show rad2deg(s.integrator[sys.elevation])


@info "System initialized at:"
toc()

# --- Stabilize system at operating point ---
@info "Stabilizing system at operating point..."
s.integrator.ps[sys.stabilize] = true
stabilization_steps = Int(10 ÷ dt)
for i in 1:stabilization_steps
    next_step!(s; dt, vsm_interval=0.05÷dt)
end
s.integrator.ps[sys.stabilize] = false

# --- Linearize at operating point ---
@info "Linearizing system at operating point..."
@time (; A, B, C, D) = KiteModels.linearize(s)
@time (; A, B, C, D) = KiteModels.linearize(s)
@show norm(A)
@info "System linearized with matrix dimensions:" A=size(A) B=size(B) C=size(C) D=size(D)

@show rad2deg(s.lin_prob[sys.elevation])

# --- Get operating point values ---
# Extract state and input at operating point
u_op = copy(s.integrator[sys.set_values])
# Create a SysState to capture state values at operating point
sys_state_op = KiteModels.SysState(s)
# Also get the direct state vector for linear model
x_op = copy(s.integrator.u)

# --- Create discretized system matrices ---
# Function to compute discretized system matrices using matrix exponential method
function discretize_linear_system(A, B, dt)
    n = size(A, 1)
    m = size(B, 2)
    
    # Create augmented matrix for simultaneous computation
    M = [A B; zeros(m, n) zeros(m, m)]
    
    # Compute matrix exponential
    expM = exp(M * dt)
    
    # Extract discretized matrices
    Ad = expM[1:n, 1:n]
    Bd = expM[1:n, n+1:n+m]
    
    return Ad, Bd
end

# Discretize the continuous-time system for more accurate integration
@info "Discretizing system matrices..."
Ad, Bd = discretize_linear_system(A, B, dt)

# Verify operating point stability
derivatives_at_op = A * zeros(size(A, 1)) + B * zeros(size(B, 2))
@info "Derivatives at operating point (should be near zero):" norm(derivatives_at_op)

# Create loggers
logger_nonlinear = Logger(length(s.point_system.points), steps)
logger_linear = deepcopy(logger_nonlinear)  # Same structure for comparison

# Initialize system state trackers
sys_state_nonlinear = KiteModels.SysState(s)
# For the linear model, we'll use a deviation from operating point
sys_state_linear = deepcopy(sys_state_op)

# --- Prepare the simulation ---
simulation_time_points = zeros(Float64, steps)
# Pre-allocate arrays with fixed size instead of growing them dynamically
input_history = Vector{Vector{Float64}}(undef, steps)
perturbation_history = Vector{Vector{Float64}}(undef, steps)

# --- Set up linear state tracking ---
# Create a linear model state vector that starts at zero (deviations from op point)
x_linear = zeros(length(x_op))
# Create output vector for linear model
y_linear = zeros(size(C, 1))

@info "Starting side-by-side simulation..."
# Begin simulation with fixed number of steps
try
    for i in 1:steps
        global x_linear, sim_time
        # Calculate current simulation time
        sim_time = (i-1) * dt
        simulation_time_points[i] = sim_time
        
        # --- Calculate time-varying steering inputs ---
        # Use sinusoidal input for more realistic testing
        steering = steering_magnitude * sin(2π * steering_freq * sim_time)
        
        # --- Nonlinear system simulation ---
        # Compute control inputs: base winch force + steering
        set_values_nonlinear = copy(u_op)
        set_values_nonlinear .+= [0.0, steering, -steering]
        input_history[i] = copy(set_values_nonlinear)
        
        # Compute perturbation from operating point
        perturbation = set_values_nonlinear - u_op
        perturbation_history[i] = copy(perturbation)
        
        # Step nonlinear simulation
        (t_new, _) = next_step!(s, set_values_nonlinear; dt, vsm_interval=vsm_interval)
        
        # Log nonlinear state
        KiteModels.update_sys_state!(sys_state_nonlinear, s)
        log!(logger_nonlinear, sys_state_nonlinear)
        
        # --- Linear system simulation ---
        # Use discretized state-space model
        # x[k+1] = Ad*x[k] + Bd*u[k]
        x_linear = Ad * x_linear + Bd * perturbation
        
        # Compute linear system output
        y_linear .= C * x_linear + D * perturbation
        sys_state_linear.turn_rates = sys_state_op.turn_rates .+ y_linear
        
        # Log linear state
        log!(logger_linear, sys_state_linear)
    end
catch e
    if isa(e, AssertionError)
        @show i
        println(e)
    else
        rethrow(e)
    end
end

@info "Simulation completed:"
toc()

# --- Save logs ---
save_log(logger_nonlinear, "nonlinear_model")
save_log(logger_linear, "linear_model")
lg_nonlinear = load_log("nonlinear_model")
lg_linear = load_log("linear_model")
sl_nonlinear = lg_nonlinear.syslog
sl_linear = lg_linear.syslog

# --- Plot comparison results ---
# Extract necessary data
turn_rates_nonlinear_deg = rad2deg.(hcat(sl_nonlinear.turn_rates...))
turn_rates_linear_deg = rad2deg.(hcat(sl_linear.turn_rates...))

# Format input data for plotting
steering_inputs = [inputs[2] - u_op[2] for inputs in input_history]

# Create comparison plots using plotx
t_plot = sl_nonlinear.time

# Prepare the data in the format expected by plotx
ω_x_comparison = [turn_rates_nonlinear_deg[1,:], turn_rates_linear_deg[1,:]]
ω_y_comparison = [turn_rates_nonlinear_deg[2,:], turn_rates_linear_deg[2,:]]
ω_z_comparison = [turn_rates_nonlinear_deg[3,:], turn_rates_linear_deg[3,:]]
steering_series = [steering_inputs]

# Create the comparison plot
p_comparison = plotx(t_plot,
    ω_x_comparison,
    ω_y_comparison,
    ω_z_comparison,
    steering_series;
    ylabels=["ω_x [°/s]", "ω_y [°/s]", "ω_z [°/s]", "Steering [Nm]"],
    labels=[
        ["Nonlinear", "Linear"],
        ["Nonlinear", "Linear"],
        ["Nonlinear", "Linear"],
        ["Input"]
    ],
    fig="Linear vs Nonlinear Model Comparison")
display(p_comparison)

# --- Calculate error metrics ---
# Function to calculate normalized RMSE
function calculate_nrmse(actual, predicted)
    # Trim to same length if needed
    min_length = min(length(actual), length(predicted))
    # Calculate RMSE
    rmse = sqrt(sum((actual[1:min_length] - predicted[1:min_length]).^2) / min_length)
    # Normalize by range of actual values
    range_actual = maximum(actual[1:min_length]) - minimum(actual[1:min_length])
    # Avoid division by zero
    if abs(range_actual) < 1e-10
        return rmse  # Return unnormalized if range is too small
    else
        return rmse / range_actual
    end
end

# Calculate both average L2 norm and NRMSE
len = turn_rates_nonlinear_deg[1,:]
error_ω_x = norm(turn_rates_nonlinear_deg[1,:] - turn_rates_linear_deg[1,:]) / len
error_ω_y = norm(turn_rates_nonlinear_deg[2,:] - turn_rates_linear_deg[2,:]) / len
error_ω_z = norm(turn_rates_nonlinear_deg[3,:] - turn_rates_linear_deg[3,:]) / len

nrmse_ω_x = calculate_nrmse(turn_rates_nonlinear_deg[1,:], turn_rates_linear_deg[1,:])
nrmse_ω_y = calculate_nrmse(turn_rates_nonlinear_deg[2,:], turn_rates_linear_deg[2,:])
nrmse_ω_z = calculate_nrmse(turn_rates_nonlinear_deg[3,:], turn_rates_linear_deg[3,:])

@info "Error metrics (average L2 norm):" error_ω_x error_ω_y error_ω_z
@info "Normalized RMSE (lower is better):" nrmse_ω_x nrmse_ω_y nrmse_ω_z

# --- Plot error over time ---
# Calculate difference between linear and nonlinear models
diff_ω_x = turn_rates_linear_deg[1,:] - turn_rates_nonlinear_deg[1,:]
diff_ω_y = turn_rates_linear_deg[2,:] - turn_rates_nonlinear_deg[2,:]
diff_ω_z = turn_rates_linear_deg[3,:] - turn_rates_nonlinear_deg[3,:]

# Plot the differences
p_error = plotx(t_plot,
    [diff_ω_x],
    [diff_ω_y],
    [diff_ω_z],
    [steering_inputs];
    ylabels=["ω_x error [°/s]", "ω_y error [°/s]", "ω_z error [°/s]", "Steering [Nm]"],
    labels=[
        ["Error"],
        ["Error"],
        ["Error"],
        ["Input"]
    ],
    fig="Linear vs Nonlinear Model Error")
display(p_error)