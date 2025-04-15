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
using ModelingToolkit: setu, getu, setp

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
set.quasi_static = true
set.bridle_fracs = [0.0, 0.93]

wing = RamAirWing(set; prn=false, n_groups=2)
aero = BodyAerodynamics([wing])
vsm_solver = Solver(aero; solver_type=NONLIN, atol=2e-8, rtol=2e-8)
point_system = create_simple_ram_point_system(set, wing)
s = RamAirKite(set, aero, vsm_solver, point_system)

measure = Measurement()
measure.set_values .= [-55, -4.0, -4.0]  # Set values of the torques of the three winches. [Nm]
set_values = measure.set_values
s.set.abs_tol = 1e-5
s.set.rel_tol = 1e-5

# Initialize at elevation
measure.sphere_pos .= deg2rad.([83.0 83.0; 1.0 -1.0])
KiteModels.init_sim!(s, measure; remake=false)
sys = s.sys
prob = s.prob
solver = FBDF(nlsolve=OrdinaryDiffEqNonlinearSolve.NLNewton(relax=0.4, max_iter=1000))
@time solve(s.prob, solver; dt, abstol=s.set.abs_tol, reltol=s.set.rel_tol, save_on=false, save_everystep=false)
@time solve(s.prob, solver; dt, abstol=s.set.abs_tol, reltol=s.set.rel_tol, save_on=false, save_everystep=false)

# Function to step simulation with input u
function step_with_input(x, u, _, p)
    (s, stiff_x, set_x, set_sx, set_u, get_x, dt) = p
    set_x(s.prob, x)
    set_sx(s.prob, stiff_x)
    set_u(s.prob, u)
    sol = solve(s.prob, solver; dt, abstol=s.set.abs_tol, reltol=s.set.rel_tol, save_on=false, save_everystep=false, save_start=false)
    return get_x(sol)[1]
end

# Get initial state
x_vec = KiteModels.get_nonstiff_unknowns(s)
sx_vec = KiteModels.get_stiff_unknowns(s)
set_x = setu(s.prob, Initial.(x_vec))
set_sx = setu(s.prob, sx_vec)
set_u = setu(s.prob, collect(sys.set_values))
get_x = getu(s.prob, x_vec)
get_sx = getu(s.prob, sx_vec)
x0 = get_x(s.prob)
sx0 = get_sx(s.prob)

# Test steering inputs and record angular velocity response
function test_response(s, input_range, input_idx; steps=1)
    angular_vels = zeros(3, length(input_range))
    total_time = 0.0
    iter = 0

    for (i, input_val) in enumerate(input_range)
        u = copy(measure.set_values)
        u[input_idx] += input_val
        x = copy(x0)
        for i in 1:steps
            p = (s, sx0, set_x, set_sx, set_u, get_x, dt)
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

