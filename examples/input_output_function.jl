#=
This example analyzes the input-output behavior of a ram air kite model by:

1. Creating a ram air kite model and initializing it at 60° elevation
2. Stabilizing the system
3. Testing the system response to different control inputs:
   - Left steering line torque
   - Right steering line torque

The script plots relationships between steering inputs and resulting angular velocities,
providing insight into the kite's steering behavior and control characteristics.
=#

using Timers
tic()
using KiteModels, LinearAlgebra, OrdinaryDiffEqCore, OrdinaryDiffEqNonlinearSolve, OrdinaryDiffEqBDF
using ModelingToolkit
using SciMLBase: successful_retcode
using ModelingToolkit: setu, getu

PLOT = true
if PLOT
    using ControlPlots
end

include(joinpath(@__DIR__, "plotting.jl"))

# Initialize model
set = load_settings("system_ram.yaml")
set.abs_tol = 0.1
set.rel_tol = 0.1
set.segments = 2
set.quasi_static = true
set.physical_model = "ram"
set.sample_freq = 20
dt = 1/set.sample_freq

s = RamAirKite(set)

measure = Measurement()
measure.set_values .= [-50, -1.0, -1.0]  # Set values of the torques of the three winches. [Nm]
set_values = measure.set_values

# Initialize at elevation
measure.sphere_pos .= deg2rad.([70.0 70.0; 1.0 -1.0])
KiteModels.init_sim!(s, measure; 
    adaptive=false, remake=false, reload=true, 
    solver=FBDF(nlsolve=OrdinaryDiffEqNonlinearSolve.NLNewton(relax=0.8, max_iter=1000))
)
OrdinaryDiffEqCore.set_proposed_dt!(s.integrator, dt)

sys = s.sys
# # Stabilize system
s.integrator.ps[sys.stabilize] = true
for i in 1:10÷dt
    next_step!(s; dt, vsm_interval=1)
end
s.integrator.ps[sys.stabilize] = false
plot(s, 0.0; zoom=true, front=false)

# Function to step simulation with input u
function step_with_input_integ(x, u, _, p)
    (s, set_x, set_ix, set_sx, sx, set_u, get_x, dt) = p
    set_x(s.integrator, x)
    set_sx(s.integrator, sx)
    set_u(s.integrator, u)
    OrdinaryDiffEqCore.set_t!(s.integrator, 0.0)
    OrdinaryDiffEqCore.reinit!(s.integrator, s.integrator.u; reinit_dae=false)
    OrdinaryDiffEqCore.step!(s.integrator, dt)
    return get_x(s.integrator)
end

# Get initial state
x_vec = KiteModels.get_nonstiff_unknowns(s)
sx_vec = KiteModels.get_stiff_unknowns(s)
set_ix = setu(s.integrator, Initial.(x_vec)) # set_ix might not be needed anymore if prob is removed
set_x = setu(s.integrator, x_vec)
set_sx = setu(s.integrator, sx_vec)
set_u = setu(s.integrator, collect(sys.set_values))
get_x = getu(s.integrator, x_vec)
get_sx = getu(s.integrator, sx_vec)
sx = get_sx(s.integrator)
x0 = get_x(s.integrator)

function test_response(s, input_range, input_idx, step_fn, u0, x_idxs=nothing; steps=1)
    # If no x_idxs specified, default to angular velocities
    if isnothing(x_idxs)
        x_idxs = (length(x0)-8):(length(x0)-6) # Assuming ω_b are the last 3 non-stiff states before stiff ones
        output_size = 3
    else
        output_size = length(x_idxs)
    end

    output = zeros(output_size, length(input_range))
    total_time = 0.0
    iter = 0

    for (i, input_val) in enumerate(input_range)
        u = copy(u0)
        u[input_idx] += input_val
        x = copy(x0)
        sx_ = copy(sx)
        set_ix(s.integrator, x)
        OrdinaryDiffEqCore.reinit!(s.integrator)
        for _ in 1:steps # Use _ if i is not used inside the loop
            p = (s, set_x, set_ix, set_sx, sx_, set_u, get_x, dt) # Pass set_ix even if unused by integ function
            total_time += @elapsed x = step_fn(x, u, nothing, p)
            iter += 1
        end
        output[:, i] = x[x_idxs]
    end

    times_rt = dt*iter/total_time
    @info "Number of steps: $iter, Times realtime: $times_rt, Total time: $total_time"
    return input_range, output, times_rt
end

# Add helper function to find state indices
function find_state_index(x_vec, symbol)
    # Compare the variables using isequal for symbolic equality
    idx = findfirst(x -> isequal(x, symbol), x_vec)
    isnothing(idx) && error("Symbol $symbol not found in state vector")
    return idx
end

function plot_input_output_relations(step_fn, suffix)
    # Find relevant state indices
    ω_idxs = [find_state_index(x_vec, sys.ω_b[i]) for i in 1:3]
    twist_idx = find_state_index(sx_vec, sys.free_twist_angle[1])

    # Test ranges
    steer_range = range(-1e-2, 1e-2, length=100)
    twist_range = range(-1e-2, 1e-2, length=100)

    # Test steering input vs omega
    @info "Testing steering input response for $suffix..."
    _, ω_steer_left, _ = test_response(s, steer_range, 2, step_fn, measure.set_values, ω_idxs)
    _, ω_steer_right, _ = test_response(s, steer_range, 3, step_fn, measure.set_values, ω_idxs)

    # Test twist angle vs omega
    @info "Testing twist angle response for $suffix..."
    function step_with_twist(x, twist_val, _, p)
        p[5][twist_idx] = twist_val[1]  # Set twist angle directly in stiff state vector copy
        return step_fn(x, measure.set_values, nothing, p)
    end
    _, ω_twist, _ = test_response(s, twist_range, 1, step_with_twist, zeros(3), ω_idxs) # u0 is zeros(3) here, seems ok

    return ω_steer_left, ω_steer_right, ω_twist, steer_range, twist_range
end

ω_steer_left_integ, ω_steer_right_integ, ω_twist_integ, steer_range, twist_range =
    plot_input_output_relations(step_with_input_integ, "integ")

# --- Plot combined results ---
# Simplified data collection
steering_data = [ω_steer_left_integ[1,:], ω_steer_right_integ[1,:]]
steering_labels = ["Left Steering (integ)", "Right Steering (integ)"]
twist_data = [ω_twist_integ[1,:]]
twist_labels = ["Twist Input (integ)"]

# Steering Plot
steering_plot = plotx(steer_range,
    steering_data, # ω_b[1] data
    [ω_steer_left_integ[2,:], ω_steer_right_integ[2,:]], # ω_b[2] data
    [ω_steer_left_integ[3,:], ω_steer_right_integ[3,:]]; # ω_b[3] data
    ylabels=["ω_b[1]", "ω_b[2]", "ω_b[3]"],
    labels=[steering_labels, steering_labels, steering_labels], # Repeat labels for each subplot
    fig="Steering Input vs Angular Velocity (Integrator)",
    xlabel="Steering Input [Nm]")

# Twist Plot
twist_plot = plotx(rad2deg.(twist_range),
    twist_data, # ω_b[1] data
    [ω_twist_integ[2,:]], # ω_b[2] data
    [ω_twist_integ[3,:]]; # ω_b[3] data
    ylabels=["ω_b[1]", "ω_b[2]", "ω_b[3]"],
    labels=[twist_labels, twist_labels, twist_labels], # Repeat labels for each subplot
    fig="Twist Angle vs Angular Velocity (Integrator)",
    xlabel="Twist Angle [deg]")

display(steering_plot)
display(twist_plot)
