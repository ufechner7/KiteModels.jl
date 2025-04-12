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
using KiteModels, LinearAlgebra, OrdinaryDiffEqCore
using ModelingToolkit
using ModelingToolkit: setu, getu

PLOT = true
if PLOT
    using ControlPlots
end

include(joinpath(@__DIR__, "plotting.jl"))

# Simulation parameters
dt = 0.05

# Initialize model
set = se("system_ram.yaml")
set.segments = 2
set_values = [-50, 0.0, 0.0]  # Initial values
set.quasi_static = true

wing = RamAirWing(set)
aero = BodyAerodynamics([wing])
vsm_solver = Solver(aero; solver_type=NONLIN, atol=1e-8, rtol=1e-8)
point_system = PointMassSystem(set, wing)
s = RamAirKite(set, aero, vsm_solver, point_system)

measure = Measurement()
s.set.abs_tol = 1e-5
s.set.rel_tol = 1e-3

# Initialize at elevation
measure.sphere_pos .= deg2rad.([60.0 60.0; 1.0 -1.0])
KiteModels.init_sim!(s, measure)
sys = s.sys

# Stabilize system
s.integrator.ps[sys.steady] = true
next_step!(s; dt=10.0, vsm_interval=1)
s.integrator.ps[sys.steady] = false

# Function to step simulation with input u
function step_with_input(x, u, _, p)
    (s, stiff_x, set_x, set_sx, get_x, dt) = p
    set_x(s.integrator, x)
    set_sx(s.integrator, stiff_x)
    OrdinaryDiffEqCore.reinit!(s.integrator, s.integrator.u; reinit_dae=false)
    next_step!(s, u; dt=dt, vsm_interval=0)
    return get_x(s.integrator)
end

# Get initial state
x_vec = KiteModels.get_nonstiff_unknowns(s)
sx_vec = KiteModels.get_stiff_unknowns(s)
set_x = setu(s.integrator, x_vec)
set_sx = setu(s.integrator, sx_vec)
get_x = getu(s.integrator, x_vec)
get_sx = getu(s.integrator, sx_vec)
x0 = get_x(s.integrator)
sx0 = get_sx(s.integrator)

# Test steering inputs and record angular velocity response
function test_response(s, input_range, input_idx; steps=1)
    angular_vels = zeros(3, length(input_range))
    total_time = 0.0
    iter = 0

    for (i, input_val) in enumerate(input_range)
        u = [-50.0, 0.0, 0.0]
        u[input_idx] = input_val
        x = copy(x0)
        for i in 1:steps
            p = (s, sx0, set_x, set_sx, get_x, dt)
            total_time += @elapsed x = step_with_input(x, u, nothing, p)
            iter += 1
        end
        angular_vels[:, i] = x[11:13]
    end
    
    times_rt = dt*iter/total_time
    @info "Number of steps: $iter, Times realtime: $times_rt"
    return input_range, angular_vels, times_rt
end

# Test left and right steering inputs
left_range = range(-1.0, 1.0, length=20)
time_vec_left, angular_vels_left, _ = test_response(s, left_range, 2)

right_range = range(-1.0, 1.0, length=20)
time_vec_right, angular_vels_right, _ = test_response(s, right_range, 3)

# Compare steering inputs effect on angular velocity
left_vs_right = plotx(time_vec_left, 
    [angular_vels_left[1,:], angular_vels_right[1,:]],
    [angular_vels_left[2,:], angular_vels_right[2,:]],
    [angular_vels_left[3,:], angular_vels_right[3,:]];
    ylabels=["ω_b[1]", "ω_b[2]", "ω_b[3]"], 
    labels=[
        ["Left Steering", "Right Steering"],
        ["Left Steering", "Right Steering"],
        ["Left Steering", "Right Steering"],
    ],
    fig="Steering Input vs Angular Velocity",
    xlabel="Steering Input Value")

display(left_vs_right)
