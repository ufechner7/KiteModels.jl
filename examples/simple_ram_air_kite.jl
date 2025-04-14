using KiteModels, LinearAlgebra

PLOT = true
if PLOT
    using Pkg
    if ! ("LaTeXStrings" ∈ keys(Pkg.project().dependencies))
        using TestEnv; TestEnv.activate()
    end
    using ControlPlots, LaTeXStrings
end

include(joinpath(@__DIR__, "plotting.jl"))

# Simulation parameters
dt = 0.05
total_time = 10  # Longer simulation to see oscillations
vsm_interval = 2
steps = Int(round(total_time / dt))

# Steering parameters
steering_freq = 1/2  # Hz - full left-right cycle frequency
steering_magnitude = 5.0      # Magnitude of steering input [Nm]

# Initialize model
set = se("system_ram.yaml")
set.segments = 3
set_values = [-50, -1.0, -1.0]  # Set values of the torques of the three winches. [Nm]
set.quasi_static = false
set.bridle_fracs = [0.088, 0.58, 0.93]
@show set.bridle_fracs

wing = RamAirWing(set; prn=false, n_groups=2)
aero = BodyAerodynamics([wing])
vsm_solver = Solver(aero; solver_type=NONLIN, atol=2e-8, rtol=2e-8)
point_system = generate_simple_ram_point_system(set, wing)
s = RamAirKite(set, aero, vsm_solver, point_system)

measure = Measurement()
s.set.abs_tol = 1e-5
s.set.rel_tol = 1e-4

# Initialize at elevation
measure.sphere_pos .= deg2rad.([60.0 60.0; 1.0 -1.0])
KiteModels.init_sim!(s, measure; remake=true)
plot(s, 0.0; zoom=true, front=false)
sys = s.sys

# # Stabilize system
# s.integrator.ps[sys.steady] = true
# next_step!(s; dt=10.0, vsm_interval=1)
# s.integrator.ps[sys.steady] = false

# logger = Logger(length(s.point_system.points), steps)
# sys_state = KiteModels.SysState(s)
# t = 0.0
# runtime = 0.0
# integ_runtime = 0.0

# try
#     while t < total_time
#         global t, set_values, runtime, integ_runtime
#         PLOT && plot(s, t; zoom=false, front=false)
        
#         # Calculate steering inputs based on cosine wave
#         steering = steering_magnitude * cos(2π * steering_freq * t+0.1)
#         set_values = -s.set.drum_radius .* s.integrator[sys.winch_force]
#         if t > 1.0
#             set_values .+= [0.0, steering, -steering]  # Opposite steering for left/right
#         end
        
#         # Step simulation
#         steptime = @elapsed (t_new, integ_steptime) = next_step!(s, set_values; dt, vsm_interval)
#         t = t_new - 10.0  # Adjust for initial stabilization time
        
#         # Track performance after initial transient
#         if (t > total_time/2)
#             runtime += steptime
#             integ_runtime += integ_steptime
#         end
        
#         # Log state variables
#         KiteModels.update_sys_state!(sys_state, s)
#         sys_state.var_01 = s.integrator[sys.ω_b[1]]
#         sys_state.var_02 = s.integrator[sys.ω_b[2]]
#         sys_state.var_03 = s.integrator[sys.ω_b[3]]

#         sys_state.var_04 = s.integrator[sys.tether_vel[2]]
#         sys_state.var_05 = s.integrator[sys.tether_vel[3]]

#         sys_state.var_06 = s.integrator[sys.aero_force_b[3]]
#         sys_state.var_07 = s.integrator[sys.aero_moment_b[2]]
#         sys_state.var_08 = s.integrator[sys.group_aero_moment[1]]

#         sys_state.var_09 = s.integrator[sys.twist_angle[1]]
#         sys_state.var_10 = s.integrator[sys.twist_angle[2]]
#         sys_state.var_11 = s.integrator[sys.twist_angle[3]]
#         sys_state.var_12 = s.integrator[sys.twist_angle[4]]

#         sys_state.var_13 = s.integrator[sys.pulley_l0[1]]
#         sys_state.var_14 = s.integrator[sys.pulley_l0[2]]
        
#         sys_state.var_15 = rad2deg(calc_aoa(s))
        
#         log!(logger, sys_state)
#     end
# catch e
#     if isa(e, AssertionError)
#         @show t
#         println(e)
#     else
#         rethrow(e)
#     end
# end

# # Plot results
# c = collect
# save_log(logger, "tmp")
# lg =load_log("tmp")
# sl = lg.syslog

# p = plotx(sl.time .- 10, 
#     [rad2deg.(sl.var_01), rad2deg.(sl.var_02), rad2deg.(sl.var_03)],
#     [c(sl.var_04), c(sl.var_05)],
#     [c(sl.var_06), c(sl.var_07), c(sl.var_08)],
#     [rad2deg.(c(sl.var_09)), rad2deg.(c(sl.var_10)), rad2deg.(c(sl.var_11)), rad2deg.(c(sl.var_12))],
#     [c(sl.var_13), c(sl.var_14)],
#     [c(sl.var_15)],
#     [rad2deg.(c(sl.heading))];
#     ylabels=["turn rates [°/s]", L"v_{ro}~[m/s]", "vsm", "twist [°]", "pulley", "AoA [°]", "heading [°]"],
#     ysize=10,
#     labels=[
#         [L"ω_x", L"ω_y", L"ω_z"],
#         ["vel[1]", "vel[2]"],
#         ["force[3]", "kite moment[2]", "group moment[1]"],
#         ["twist_angle[1]", "twist_angle[2]", "twist_angle[3]", "twist_angle[4]"],
#         ["pulley_l0[1]", "pulley_l0[2]"],
#         ["angle of attack"],
#         ["heading"]
#     ],
#     fig="Oscillating Steering Input Response")
# display(p)

# @info "Performance:" times_realtime=(total_time/2)/runtime integrator_times_realtime=(total_time/2)/integ_runtime
