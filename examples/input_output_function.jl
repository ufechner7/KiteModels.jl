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
set.abs_tol = 5e-5
set.rel_tol = 5e-5
set.segments = 2
set.quasi_static = false
set.physical_model = "simple_ram"
if set.physical_model == "ram"
    set.bridle_fracs = [0.088, 0.31, 0.58, 0.93]
elseif set.physical_model == "simple_ram"
    set.bridle_fracs = [0.0, 0.93]
end
set.sample_freq = 20
dt = 1/set.sample_freq

s = RamAirKite(set)

measure = Measurement()
measure.set_values .= [-50, -1.0, -1.0]  # Set values of the torques of the three winches. [Nm]
set_values = measure.set_values

# Initialize at elevation
measure.sphere_pos .= deg2rad.([60.0 60.0; 1.0 -1.0])
KiteModels.init_sim!(s, measure; remake=false, reload=true)

sys = s.sys

# Function to step simulation with input u
function step_with_input_integ(x, u, _, p)
    (s, set_x, set_ix, set_sx, sx, set_u, get_x, dt) = p
    set_x(s.integrator, x)
    set_sx(s.integrator, sx)
    set_u(s.integrator, u)
    OrdinaryDiffEqCore.reinit!(s.integrator, s.integrator.u; reinit_dae=false)
    OrdinaryDiffEqCore.step!(s.integrator, dt)
    return get_x(s.integrator)
end

solver = FBDF(nlsolve=OrdinaryDiffEqNonlinearSolve.NLNewton(relax=0.4))

# Function to step simulation with input u
function step_with_input_prob(x, u, _, p)
    (s, set_x, set_ix, set_sx, sx, set_u, get_x, dt) = p
    set_ix(s.prob, x)
    set_u(s.prob, u)
    sol = solve(s.prob, solver; dt, abstol=s.set.abs_tol, reltol=s.set.rel_tol, save_on=false, save_everystep=false, save_start=false)
    return get_x(sol)[1]
end

measure.sphere_pos .= deg2rad.([83.0 83.0; 1.0 -1.0])
KiteModels.reinit!(s, measure; reload=false)

# Get initial state
x_vec = KiteModels.get_nonstiff_unknowns(s)
sx_vec = KiteModels.get_stiff_unknowns(s)
set_ix = setu(s.integrator, Initial.(x_vec))
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
        x_idxs = (length(x0)-8):(length(x0)-6)
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
        for i in 1:steps
            p = (s, set_x, set_ix, set_sx, sx_, set_u, get_x, dt)
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
    twist_idx = find_state_index(x_vec, sys.free_twist_angle[1])
    
    # Test ranges
    steer_range = range(-0.1, 0.1, length=100)
    twist_range = range(-0.1, 0.1, length=100)
    
    # Test steering input vs omega
    @info "Testing steering input response for $suffix..."
    _, ω_steer_left, _ = test_response(s, steer_range, 2, step_fn, measure.set_values, ω_idxs)
    _, ω_steer_right, _ = test_response(s, steer_range, 3, step_fn, measure.set_values, ω_idxs)

    # Test twist angle vs omega 
    @info "Testing twist angle response for $suffix..."
    function step_with_twist(x, twist_val, _, p)
        x[twist_idx] = twist_val[1]  # Set twist angle directly
        return step_fn(x, measure.set_values, nothing, p)
    end
    _, ω_twist, _ = test_response(s, twist_range, 1, step_with_twist, zeros(3), ω_idxs)

    return ω_steer_left, ω_steer_right, ω_twist, steer_range, twist_range
end

# Run analysis for both methods
ω_steer_left_prob, ω_steer_right_prob, ω_twist_prob, steer_range, twist_range = 
    plot_input_output_relations(step_with_input_prob, "prob")
ω_steer_left_integ, ω_steer_right_integ, ω_twist_integ, _, _ = 
    plot_input_output_relations(step_with_input_integ, "integ")

# Plot combined results
steering_plot = plotx(steer_range, 
    [ω_steer_left_prob[1,:], ω_steer_right_prob[1,:], 
     ω_steer_left_integ[1,:], ω_steer_right_integ[1,:]],
    [ω_steer_left_prob[2,:], ω_steer_right_prob[2,:],
     ω_steer_left_integ[2,:], ω_steer_right_integ[2,:]],
    [ω_steer_left_prob[3,:], ω_steer_right_prob[3,:],
     ω_steer_left_integ[3,:], ω_steer_right_integ[3,:]];
    ylabels=["ω_b[1]", "ω_b[2]", "ω_b[3]"], 
    labels=[
        ["Left Steering (prob)", "Right Steering (prob)", 
         "Left Steering (integ)", "Right Steering (integ)"],
        ["Left Steering (prob)", "Right Steering (prob)",
         "Left Steering (integ)", "Right Steering (integ)"],
        ["Left Steering (prob)", "Right Steering (prob)",
         "Left Steering (integ)", "Right Steering (integ)"],
    ],
    fig="Steering Input vs Angular Velocity Comparison",
    xlabel="Steering Input [Nm]")

twist_plot = plotx(rad2deg.(twist_range),
    [ω_twist_prob[1,:], ω_twist_integ[1,:]],
    [ω_twist_prob[2,:], ω_twist_integ[2,:]],
    [ω_twist_prob[3,:], ω_twist_integ[3,:]];
    ylabels=["ω_b[1]", "ω_b[2]", "ω_b[3]"],
    labels=[
        ["Twist Input (prob)", "Twist Input (integ)"],
        ["Twist Input (prob)", "Twist Input (integ)"],
        ["Twist Input (prob)", "Twist Input (integ)"]
    ],
    fig="Twist Angle vs Angular Velocity Comparison",
    xlabel="Twist Angle [deg]")

display(steering_plot)
display(twist_plot)
