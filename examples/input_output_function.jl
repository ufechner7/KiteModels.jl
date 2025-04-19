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
using ModelingToolkit: setu, getu

PLOT = true
if PLOT
    using ControlPlots
end

include(joinpath(@__DIR__, "plotting.jl"))

# Initialize model
set = se("system_ram.yaml")
set.segments = 2
set.quasi_static = true
set.physical_model = "ram"
if set.physical_model == "ram"
    set.bridle_fracs = [0.088, 0.31, 0.58, 0.93]
elseif set.physical_model == "simple_ram"
    set.bridle_fracs = [0.0, 0.93]
end
set.sample_freq = 20
dt = 1/set.sample_freq

s = RamAirKite(set)

measure = Measurement()
measure.set_values .= [-55, -4.0, -4.0]  # Set values of the torques of the three winches. [Nm]
set_values = measure.set_values
s.set.abs_tol = 0.01
s.set.rel_tol = 0.01

# Initialize at elevation
measure.sphere_pos .= deg2rad.([83.0 83.0; 1.0 -1.0])
KiteModels.init_sim!(s, measure; remake=false)
sys = s.sys

# Stabilize system
s.integrator.ps[sys.steady] = true
next_step!(s; dt=10.0, vsm_interval=1)
s.integrator.ps[sys.steady] = false

# Function to step simulation with input u
function step_with_input_integ(x, u, _, p)
    (s, set_x, set_u, get_x, dt) = p
    set_x(s.integrator, x)
    set_u(s.integrator, u)
    OrdinaryDiffEqCore.reinit!(s.integrator, s.integrator.u; reinit_dae=true)
    OrdinaryDiffEqCore.step!(s.integrator, dt)
    return get_x(s.integrator)
end

solver = FBDF(nlsolve=OrdinaryDiffEqNonlinearSolve.NLNewton(relax=0.4))

# Function to step simulation with input u
function step_with_input_prob(x, u, _, p)
    (s, set_x, set_u, get_x, dt) = p
    set_x(s.prob, x)
    set_u(s.prob, u)
    sol = solve(s.prob, solver; dt, abstol=s.set.abs_tol, reltol=s.set.rel_tol, save_on=false, save_everystep=false, save_start=false)
    return get_x(sol)[1]
end

# Get initial state
x_vec = KiteModels.get_unknowns(s)
set_x = setu(s.integrator, Initial.(x_vec))
set_u = setu(s.integrator, collect(sys.set_values))
get_x = getu(s.integrator, x_vec)
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
        for i in 1:steps
            p = (s, set_x, set_u, get_x, dt)
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

function plot_input_output_relations(step_fn)
    # Find relevant state indices
    ω_idxs = [find_state_index(x_vec, sys.ω_b[i]) for i in 1:3]
    twist_idx = find_state_index(x_vec, sys.free_twist_angle[1])
    
    # Test ranges
    steer_range = range(-0.1, 0.1, length=20)
    twist_range = range(-0.1, 0.1, length=20)
    
    # Test steering input vs omega
    @info "Testing steering input response..."
    _, ω_steer_left, _ = test_response(s, steer_range, 2, step_fn, measure.set_values, ω_idxs)
    _, ω_steer_right, _ = test_response(s, steer_range, 3, step_fn, measure.set_values, ω_idxs)

    # Test twist angle vs omega 
    @info "Testing twist angle response..."
    function step_with_twist(x, twist_val, _, p)
        x[twist_idx] = twist_val[1]  # Set twist angle directly
        return step_fn(x, measure.set_values, nothing, p)
    end
    _, ω_twist, _ = test_response(s, twist_range, 1, step_with_twist, zeros(3), ω_idxs)

    # Plot results
    steering_plot = plotx(steer_range, 
        [ω_steer_left[1,:], ω_steer_right[1,:]],
        [ω_steer_left[2,:], ω_steer_right[2,:]],
        [ω_steer_left[3,:], ω_steer_right[3,:]];
        ylabels=["ω_b[1]", "ω_b[2]", "ω_b[3]"], 
        labels=[
            ["Left Steering", "Right Steering"],
            ["Left Steering", "Right Steering"],
            ["Left Steering", "Right Steering"],
        ],
        fig="Steering Input vs Angular Velocity",
        xlabel="Steering Input [Nm]")

    twist_plot = plotx(rad2deg.(twist_range),
        [ω_twist[1,:]], [ω_twist[2,:]], [ω_twist[3,:]];
        ylabels=["ω_b[1]", "ω_b[2]", "ω_b[3]"],
        labels=[["Twist Input"], ["Twist Input"], ["Twist Input"]],
        fig="Twist Angle vs Angular Velocity",
        xlabel="Twist Angle [deg]")

    return steering_plot, twist_plot
end

# Run analysis and display plots
steer_plot, twist_plot = plot_input_output_relations(step_with_input_prob)
# steer_plot, twist_plot = plot_input_output_relations(step_with_input_integ)
display(steer_plot)
display(twist_plot)
