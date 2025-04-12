using KiteModels, LinearAlgebra

PLOT = true
if PLOT
    using ControlPlots
end

include(joinpath(@__DIR__, "plotting.jl"))

# Simulation parameters
dt = 0.05
total_time = 3.8  # Longer simulation to see oscillations
vsm_interval = 5
steps = Int(round(total_time / dt))

# Steering parameters
steering_freq = 1/2  # Hz - full left-right cycle frequency
steering_magnitude = 5.0  # Magnitude of steering input

# Initialize model
set = se("system_ram.yaml")
set.segments = 6
set_values = [-50, 0.0, 0.0]  # Initial values
set.quasi_static = false

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

logger = Logger(length(s.point_system.points), steps)
sys_state = KiteModels.SysState(s)
t = 0.0
runtime = 0.0
integ_runtime = 0.0

try
    while t < total_time
        global t, set_values, runtime, integ_runtime
        PLOT && plot(s, t; zoom=false, front=false)
        
        # Calculate steering inputs based on cosine wave
        steering = steering_magnitude * cos(2π * steering_freq * t)
        set_values = -s.set.drum_radius .* s.integrator[sys.winch_force]
        if t > 1.0
            set_values .+= [0.0, steering, -steering]  # Opposite steering for left/right
        end
        
        # Step simulation
        steptime = @elapsed (t_new, integ_steptime) = next_step!(s, set_values; dt, vsm_interval)
        t = t_new - 10.0  # Adjust for initial stabilization time
        
        # Track performance after initial transient
        if (t > total_time/2)
            runtime += steptime
            integ_runtime += integ_steptime
        end
        
        # Log state variables
        KiteModels.update_sys_state!(sys_state, s)
        sys_state.var_01 = s.integrator[sys.ω_b[1]]
        sys_state.var_02 = s.integrator[sys.ω_b[2]]
        sys_state.var_03 = s.integrator[sys.ω_b[3]]

        sys_state.var_04 = s.integrator[sys.tether_vel[2]]
        sys_state.var_05 = s.integrator[sys.tether_vel[3]]

        sys_state.var_06 = s.integrator[sys.aero_force_b[3]]
        sys_state.var_07 = s.integrator[sys.aero_moment_b[2]]
        sys_state.var_08 = s.integrator[sys.group_aero_moment[1]]

        sys_state.var_09 = s.integrator[sys.twist_angle[1]]
        sys_state.var_10 = s.integrator[sys.twist_angle[2]]
        sys_state.var_11 = s.integrator[sys.twist_angle[3]]
        sys_state.var_12 = s.integrator[sys.twist_angle[4]]

        sys_state.var_13 = clamp(s.integrator[sys.pulley_acc[1]], -100, 100)
        sys_state.var_14 = clamp(s.integrator[sys.pulley_acc[2]], -100, 100)
        
        # Calculate apparent wind angle
        va_kite_b = s.integrator[sys.va_kite_b]
        e_x = s.integrator[sys.e_x]
        sys_state.var_15 = rad2deg(acos(dot(normalize(va_kite_b), e_x)))
        
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

# Plot results
p = plotx(logger.time_vec, 
    [logger.var_01_vec, logger.var_02_vec, logger.var_03_vec],
    [logger.var_04_vec, logger.var_05_vec],
    [logger.var_06_vec, logger.var_07_vec, logger.var_08_vec],
    [logger.var_09_vec, logger.var_10_vec, logger.var_11_vec, logger.var_12_vec],
    [logger.var_13_vec, logger.var_14_vec],
    [logger.var_15_vec],
    [logger.heading_vec];
    ylabels=["kite", "tether", "vsm", "twist", "pulley", "wind", "heading"],
    labels=[
        ["ω_b[1]", "ω_b[2]", "ω_b[3]"],
        ["vel[1]", "vel[2]"],
        ["force[3]", "kite moment[2]", "group moment[1]"],
        ["twist_angle[1]", "twist_angle[2]", "twist_angle[3]", "twist_angle[4]"],
        ["pulley acc[1]", "pulley acc[2]"],
        ["angle of attack"],
        ["heading"]
    ],
    fig="Oscillating Steering Input Response")
display(p)

@info "Performance:" times_realtime=(total_time/2)/runtime integrator_times_realtime=(total_time/2)/integ_runtime
