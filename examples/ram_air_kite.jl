using Timers
tic()
@info "Loading packages "

using KiteModels, LinearAlgebra, Statistics

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

# Simulation parameters
dt = 0.05
total_time = 10  # Longer simulation to see oscillations
vsm_interval = 3
steps = Int(round(total_time / dt))

# Steering parameters
steering_freq = 1/2  # Hz - full left-right cycle frequency
steering_magnitude = 5.0      # Magnitude of steering input [Nm]

# Initialize model
set = load_settings("system_ram.yaml")
set.segments = 3
set_values = [-50, 0.0, 0.0]  # Set values of the torques of the three winches. [Nm]
set.quasi_static = false
set.physical_model = "ram"
if set.physical_model == "ram"
    set.bridle_fracs = [0.088, 0.31, 0.58, 0.93]
elseif set.physical_model == "simple_ram"
    set.l_tether = 53
    set.bridle_fracs = [0.25, 0.93]
end

@info "Creating wing, aero, vsm_solver, point_system and s:"
s = RamAirKite(set)
toc()

# init_Q_b_w, R_b_w = KiteModels.measure_to_q(measure)
# init_kite_pos = init!(s.point_system, s.set, R_b_w)
# plot(s.point_system, 0.0; zoom=false, front=true)

measure = Measurement()
s.set.abs_tol = 1e-2
s.set.rel_tol = 1e-2

# Initialize at elevation
measure.sphere_pos .= deg2rad.([80.0 80.0; 1.0 -1.0])
KiteModels.init_sim!(s, measure; remake=false, reload=true)
sys = s.sys

@info "System initialized at:"
toc()

# Stabilize system
s.integrator.ps[sys.stiffness_frac] = 0.01
s.integrator.ps[sys.steady] = true
next_step!(s; dt=10.0, vsm_interval=1)
s.integrator.ps[sys.steady] = false
s.integrator.ps[sys.ω_damp] = 1000

"""
TODO:
- initialize with azimuth and elevation, not xyz kite pos
- calculate effective power chord frac and steering line chord frac
    - balance the front 3 lines, keep the back line at the same position
    - (F_1 * x_1 + F_2 * x_2 + F_3 * x_3) / F_total
"""

logger = Logger(length(s.point_system.points), steps)
sys_state = KiteModels.SysState(s)
t = 0.0
runtime = 0.0
integ_runtime = 0.0

try
    while t < total_time
        local steering
        global t, set_values, runtime, integ_runtime
        # PLOT && plot(s, t; zoom=false, front=true)
        
        # Calculate steering inputs based on cosine wave
        steering = steering_magnitude * cos(2π * steering_freq * t+0.1)
        set_values = -s.set.drum_radius .* s.integrator[sys.winch_force]
        _vsm_interval = 1
        if t > 1.0
            set_values .+= [0.0, steering, -steering]  # Opposite steering for left/right
            _vsm_interval = vsm_interval
        end
        
        # Step simulation
        steptime = @elapsed (t_new, integ_steptime) = next_step!(s, set_values; dt, vsm_interval=_vsm_interval)
        t = t_new - 10.0  # Adjust for initial stabilization time
        
        # Track performance after initial transient
        if (t > total_time/2)
            @show dt/integ_steptime
            runtime += steptime
            integ_runtime += integ_steptime
        end
        
        # Log state variables
        KiteModels.update_sys_state!(sys_state, s)
        sys_state.var_01 = s.integrator[sys.ω_b[1]]
        sys_state.var_02 = s.integrator[sys.ω_b[2]]
        sys_state.var_03 = s.integrator[sys.ω_b[3]]

        sys_state.var_04 = s.integrator[sys.tether_length[2]]
        sys_state.var_05 = s.integrator[sys.tether_length[3]]

        sys_state.var_06 = s.integrator[sys.aero_force_b[3]]
        sys_state.var_07 = s.integrator[sys.aero_moment_b[2]]
        sys_state.var_08 = s.integrator[sys.group_aero_moment[1]]

        sys_state.var_09 = s.integrator[sys.twist_angle[1]]
        sys_state.var_10 = s.integrator[sys.twist_angle[2]]
        
        sys_state.var_11 = norm(s.integrator[sys.kite_pos])

        sys_state.var_12 = s.integrator[sys.angle_of_attack]
        
        # @show mean(s.integrator[sys.front_frac])
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
@info "Total time without plotting:"
toc()

# Plot results
c = collect
save_log(logger, "tmp")
lg =load_log("tmp")
sl = lg.syslog

if PLOT
    p = plotx(sl.time .- 10, 
        [rad2deg.(sl.var_01), rad2deg.(sl.var_02), rad2deg.(sl.var_03)],
        [c(sl.var_04), c(sl.var_05)],
        [c(sl.var_06), c(sl.var_07), c(sl.var_08)],
        [rad2deg.(c(sl.var_09)), rad2deg.(c(sl.var_10))],
        [c(sl.var_11)],
        [rad2deg.(c(sl.var_12))],
        [rad2deg.(c(sl.heading))];
        ylabels=["turn rates [°/s]", "tether length [m]", "vsm", "twist [°]", "distance [m]", "AoA [°]", "heading [°]"],
        ysize=10,
        labels=[
            [L"ω_x", L"ω_y", L"ω_z"],
            ["len[1]", "len[2]"],
            ["force[3]", "kite moment[2]", "group moment[1]"],
            ["twist_angle[1]", "twist_angle[2]"],
            ["distance"],
            ["angle of attack"],
            ["heading"]
        ],
        fig="Oscillating Steering Input Response")
    display(p)
end

@info "Performance:" times_realtime=(total_time/2)/runtime integrator_times_realtime=(total_time/2)/integ_runtime
