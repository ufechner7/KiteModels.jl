using Timers
tic()
@info "Loading packages "

using KiteModels, LinearAlgebra

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
vsm_interval = 2
steps = Int(round(total_time / dt))

# Steering parameters
steering_freq = 1/2  # Hz - full left-right cycle frequency
steering_magnitude = 5.0      # Magnitude of steering input [Nm]

# Initialize model
set = se("system_ram.yaml")
set.segments = 3
set_values = [-50, 0.0, 0.0]  # Set values of the torques of the three winches. [Nm]
set.quasi_static = false

@info "Creating wing, aero, vsm_solver, point_system and s:"
wing = RamAirWing(set; prn=false)
aero = BodyAerodynamics([wing])
vsm_solver = Solver(aero; solver_type=NONLIN, atol=2e-8, rtol=2e-8)
point_system = PointMassSystem(set, wing)
s = RamAirKite(set, aero, vsm_solver, point_system)
toc()

measure = Measurement()
s.set.abs_tol = 1e-5
s.set.rel_tol = 1e-4

# Initialize at elevation
measure.sphere_pos .= deg2rad.([60.0 60.0; 1.0 -1.0])
KiteModels.init_sim!(s, measure)
sys = s.sys

@info "System initialized at:"
toc()

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
        local steering
        global t, set_values, runtime, integ_runtime
        PLOT && plot(s, t; zoom=false, front=false)

        # Calculate steering inputs based on cosine wave
        steering = steering_magnitude * cos(2π * steering_freq * t+0.1)
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
        @show sys_state.v_reelout
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

# --- Updated Plotting ---
# Extract necessary data using meaningful names
turn_rates_deg = rad2deg.(hcat(sl.turn_rates...))
v_reelout_23 = [sl.v_reelout[i][2] for i in eachindex(sl.v_reelout)], [sl.v_reelout[i][3] for i in eachindex(sl.v_reelout)] # Winch 2 and 3
aero_force_z = [sl.aero_force_b[i][3] for i in eachindex(sl.aero_force_b)]
aero_moment_y = [sl.aero_moment_b[i][2] for i in eachindex(sl.aero_moment_b)]
twist_angles_deg = rad2deg.(hcat(sl.twist_angles...))
AoA_deg = rad2deg.(sl.AoA)
heading_deg = rad2deg.(sl.heading)

p = plotx(sl.time .- 10,
    [turn_rates_deg[1,:], turn_rates_deg[2,:], turn_rates_deg[3,:]],
    v_reelout_23,
    [aero_force_z, aero_moment_y],
    [twist_angles_deg[1,:], twist_angles_deg[2,:], twist_angles_deg[3,:], twist_angles_deg[4,:]],
    [AoA_deg],
    [heading_deg];
    ylabels=["turn rates [°/s]", L"v_{ro}~[m/s]", "aero F/M", "twist [°]", "AoA [°]", "heading [°]"],
    ysize=10,
    labels=[
        [L"\omega_x", L"\omega_y", L"\omega_z"],
        ["v_ro[2]", "v_ro[3]"],
        [L"F_{aero,z}", L"M_{aero,y}"],
        ["twist[1]", "twist[2]", "twist[3]", "twist[4]"],
        ["AoA"],
        ["heading"]
    ],
    fig="Oscillating Steering Input Response")
display(p)

@info "Performance:" times_realtime=(total_time/2)/runtime integrator_times_realtime=(total_time/2)/integ_runtime
