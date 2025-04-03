using Revise, KiteModels, LinearAlgebra, VortexStepMethod
using ModelingToolkit
using ModelingToolkit: setu, getu
using ControlPlots

dt = 0.05
total_time = 1.5
steps = Int(round(total_time / dt))

set = se("system_ram.yaml")
set.segments = 2

if !@isdefined s
    wing = RamAirWing("data/ram_air_kite_body.obj", "data/ram_air_kite_foil.dat"; mass=set.mass, crease_frac=0.82, align_to_principal=true)
    aero = BodyAerodynamics([wing])
    vsm_solver = Solver(aero; solver_type=NONLIN, atol=1e-8, rtol=1e-8)
    point_system = PointMassSystem(set, wing)
    s = RamAirKite(set, aero, vsm_solver, point_system)

    measure = Measurement()
    # KiteModels.init_sim!(s, measure)
end
s.set.abs_tol = 1e-5
s.set.rel_tol = 1e-3
measure.sphere_pos .= deg2rad.([60.0 60.0; 1.0 -1.0])


@time KiteModels.reinit!(s, measure)
sys = s.sys
s.integrator.ps[sys.steady] = true
next_step!(s; dt=10.0, vsm_interval=1)
s.integrator.ps[sys.steady] = false

x_vec = KiteModels.get_nonstiff_unknowns(s)
sx_vec = KiteModels.get_stiff_unknowns(s)
set_x = setu(s.integrator, x_vec)
set_sx = setu(s.integrator, sx_vec)
get_x = getu(s.integrator, x_vec)
get_sx = getu(s.integrator, sx_vec)

# Function to step simulation with input u
function step_with_input(x, u, _, p)
    (s, stiff_x, set_x, set_sx, get_x, dt) = p
    set_x(s.integrator, x)
    set_sx(s.integrator, stiff_x)
    @time next_step!(s; set_values=u, dt=dt, vsm_interval=0)
    return get_x(s.integrator)
end

# Get initial state
x0 = get_x(s.integrator)
sx0 = get_sx(s.integrator)

# Test a range of inputs and record responses
function test_response(s, input_range, input_idx, num_steps=3)
    # Initialize storage for results
    time_vec = zeros(length(input_range))
    
    # Initialize response matrices for each group to plot
    tether_lengths = zeros(3, length(input_range))
    tether_vels = zeros(3, length(input_range))
    orientations = zeros(4, length(input_range))
    angular_vels = zeros(3, length(input_range))
    kite_positions = zeros(3, length(input_range))
    kite_vels = zeros(3, length(input_range))

    # For each input value
    for (i, input_val) in enumerate(input_range)
        u = [-50.0, 0.0, 0.0]
        u[input_idx] = input_val
        
        # Reset to initial state
        x = copy(x0)
        
        # Run multiple steps if requested
        for step in 1:num_steps
            p = (s, sx0, set_x, set_sx, get_x, dt)
            x = step_with_input(x, u, nothing, p)
        end
        
        # Record responses
        time_vec[i] = input_val  # Use input value as "time" for x-axis
        
        # Extract state components
        tether_lengths[:, i] = [x[1], x[3], x[5]]
        tether_vels[:, i] = [x[2], x[4], x[6]]
        orientations[:, i] = x[7:10]
        angular_vels[:, i] = x[11:13]
        kite_positions[:, i] = x[14:16]
        kite_vels[:, i] = x[17:19]
    end
    
    return time_vec, tether_lengths, tether_vels, orientations, angular_vels, kite_positions, kite_vels
end

# Test power input response (first input)
power_range = range(-50.0, -40.0, length=20)
time_vec_power, tether_lengths_power, tether_vels_power, orientations_power, 
    angular_vels_power, kite_positions_power, kite_vels_power = test_response(s, power_range, 1)

# Plot power response
p_power = plotx(time_vec_power, 
    [tether_vels_power[1,:], tether_vels_power[2,:], tether_vels_power[3,:]],
    [angular_vels_power[1,:], angular_vels_power[2,:], angular_vels_power[3,:]],
    [kite_positions_power[1,:], kite_positions_power[2,:], kite_positions_power[3,:]],
    [kite_vels_power[1,:], kite_vels_power[2,:], kite_vels_power[3,:]];
    ylabels=["Tether Velocities", "Angular Velocities", "Kite Position", "Kite Velocity"], 
    labels=[
        ["Tether Vel 1", "Tether Vel 2", "Tether Vel 3"],
        ["ω_b[1]", "ω_b[2]", "ω_b[3]"],
        ["Pos X", "Pos Y", "Pos Z"],
        ["Vel X", "Vel Y", "Vel Z"],
    ],
    fig="Power Input Response",
    xlabel="Power Input (Negative = Reel Out)")

display(p_power)

# Test left steering input response (second input)
left_range = range(-10.0, 10.0, length=20)
time_vec_left, tether_lengths_left, tether_vels_left, orientations_left, 
angular_vels_left, kite_positions_left, kite_vels_left = test_response(s, left_range, 2)

# Plot left steering response
p_left = plotx(time_vec_left, 
    [tether_vels_left[1,:], tether_vels_left[2,:], tether_vels_left[3,:]],
    [angular_vels_left[1,:], angular_vels_left[2,:], angular_vels_left[3,:]],
    [orientations_left[1,:], orientations_left[2,:], orientations_left[3,:], orientations_left[4,:]],
    [kite_vels_left[1,:], kite_vels_left[2,:], kite_vels_left[3,:]];
    ylabels=["Tether Velocities", "Angular Velocities", "Orientation Quaternion", "Kite Velocity"], 
    labels=[
        ["Tether Vel 1", "Tether Vel 2", "Tether Vel 3"],
        ["ω_b[1]", "ω_b[2]", "ω_b[3]"],
        ["q[1]", "q[2]", "q[3]", "q[4]"],
        ["Vel X", "Vel Y", "Vel Z"],
    ],
    fig="Left Steering Input Response",
    xlabel="Left Steering Input")

display(p_left)

# Test right steering input response (third input)
right_range = range(-5.0, 5.0, length=20)
time_vec_right, tether_lengths_right, tether_vels_right, orientations_right, 
angular_vels_right, kite_positions_right, kite_vels_right = test_response(s, right_range, 3)

# Plot right steering response
p_right = plotx(time_vec_right, 
    [tether_vels_right[1,:], tether_vels_right[2,:], tether_vels_right[3,:]],
    [angular_vels_right[1,:], angular_vels_right[2,:], angular_vels_right[3,:]],
    [orientations_right[1,:], orientations_right[2,:], orientations_right[3,:], orientations_right[4,:]],
    [kite_vels_right[1,:], kite_vels_right[2,:], kite_vels_right[3,:]];
    ylabels=["Tether Velocities", "Angular Velocities", "Orientation Quaternion", "Kite Velocity"], 
    labels=[
        ["Tether Vel 1", "Tether Vel 2", "Tether Vel 3"],
        ["ω_b[1]", "ω_b[2]", "ω_b[3]"],
        ["q[1]", "q[2]", "q[3]", "q[4]"],
        ["Vel X", "Vel Y", "Vel Z"],
    ],
    fig="Right Steering Input Response",
    xlabel="Right Steering Input")

display(p_right)

# Compare both steering inputs effect on angular velocity
p_comparison = plotx(time_vec_left, 
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

display(p_comparison)

nothing