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

# TODO: USE SPARSE AUTODIFF, AND LINEARIZE AROUND CORRECT OPERATING POINT

include(joinpath(@__DIR__, "plotting.jl"))

# Simulation parameters
dt = 0.001
total_time = 0.1  # Seconds of simulation after stabilization
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
measure.sphere_pos .= deg2rad.([70.0 70.0; 1.0 -1.0])
KiteModels.init_sim!(s, measure; 
    remake=false,
    reload=true,
    lin_outputs=[ω_b...]  # Specify which outputs to track in linear model
)
sys = s.sys

@info "System initialized at:"
toc()

# --- Stabilize system at operating point ---
@info "Stabilizing system at operating point..."
s.integrator.ps[sys.stabilize] = true
# stabilization_steps = Int(10 ÷ dt)  # 10 seconds of stabilization
stabilization_steps = 100
for i in 1:stabilization_steps
    next_step!(s; dt, vsm_interval=1)
end
s.integrator.ps[sys.stabilize] = false

# --- Linearize at operating point ---
@info "Linearizing system at operating point..."
@time (; A, B, C, D) = KiteModels.linearize(s)
@time (; A, B, C, D) = KiteModels.linearize(s)
@show norm(A) # 2.1875118055983325e6
@info "System linearized with matrix dimensions:" A=size(A) B=size(B) C=size(C) D=size(D)

# --- Get operating point values ---
# Extract state and input at operating point
u_op = copy(s.integrator[sys.set_values])
# Create a SysState to capture state values at operating point
sys_state_op = KiteModels.SysState(s)
# Also get the direct state vector for linear model
x_op = copy(s.integrator.u)

# Create loggers
logger_nonlinear = Logger(length(s.point_system.points), steps)
logger_linear = deepcopy(logger_nonlinear)  # Same structure for comparison

# Initialize system state trackers
sys_state_nonlinear = KiteModels.SysState(s)
# For the linear model, we'll use a deviation from operating point
sys_state_linear = deepcopy(sys_state_op)

# --- Prepare the simulation ---
sim_time = 0.0
simulation_time_points = Float64[]
# Input history for plotting
input_history = Vector{Float64}[]
# Perturbation input history
perturbation_history = Vector{Float64}[]

# --- Set up linear state tracking ---
# Create a linear model state vector that starts at zero (deviations from op point)
x_linear = zeros(length(x_op))
# Create output vector for linear model
y_linear = zeros(size(C, 1))
# Linear system: dx/dt = A*x + B*u
# We'll use a simple Euler method to integrate it alongside the nonlinear system

@info "Starting side-by-side simulation..."
# Begin simulation
try
    global sim_time
    sim_time = 0.0  # Start at t=0 for the comparison

    while sim_time < total_time
        push!(simulation_time_points, sim_time)
        
        # --- Calculate steering inputs ---
        steering = 10.0
        
        # --- Nonlinear system simulation ---
        # Compute control inputs: base winch force + steering
        set_values_nonlinear = copy(u_op)
        set_values_nonlinear .+= [0.0, steering, -steering]
        push!(input_history, copy(set_values_nonlinear))
        
        # Compute perturbation from operating point
        perturbation = set_values_nonlinear - u_op
        push!(perturbation_history, copy(perturbation))
        
        # Step nonlinear simulation
        (t_new, _) = next_step!(s, set_values_nonlinear; dt, vsm_interval=vsm_interval)
        
        # Update time
        sim_time = t_new - stabilization_steps*dt  # Adjust for initial stabilization time
        
        # Log nonlinear state
        KiteModels.update_sys_state!(sys_state_nonlinear, s)
        log!(logger_nonlinear, sys_state_nonlinear)
        
        # --- Linear system simulation ---
        # Linear system step using Euler integration
        # Perturbation input: du = u - u_op
        # dx/dt = A*x + B*du
        # For dt small enough, x(t+dt) ≈ x(t) + dx/dt * dt
        dx_dt = A * x_linear + B * perturbation
        x_linear .+= dx_dt * dt
        
        # Compute linear system output
        y_linear .= C * x_linear + D * perturbation
        sys_state_linear.turn_rates = sys_state_op.turn_rates .+ y_linear
        
        # Log linear state
        log!(logger_linear, sys_state_linear)
    end
catch e
    if isa(e, AssertionError)
        @show sim_time
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
# For each turn rate component, we'll have [nonlinear_series, linear_series]
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
# Trim to same length if needed
min_length = min(length(turn_rates_nonlinear_deg[1,:]), length(turn_rates_linear_deg[1,:]))
error_ω_x = norm(turn_rates_nonlinear_deg[1,1:min_length] - turn_rates_linear_deg[1,1:min_length]) / min_length
error_ω_y = norm(turn_rates_nonlinear_deg[2,1:min_length] - turn_rates_linear_deg[2,1:min_length]) / min_length
error_ω_z = norm(turn_rates_nonlinear_deg[3,1:min_length] - turn_rates_linear_deg[3,1:min_length]) / min_length

@info "Error metrics (average L2 norm):" error_ω_x error_ω_y error_ω_z

