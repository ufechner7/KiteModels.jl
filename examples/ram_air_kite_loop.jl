using KiteModels, LinearAlgebra

PLOT = true
if PLOT
    using ControlPlots
end

include("./plotting.jl")

dt = 0.05
total_time = 4.5
vsm_interval = 5
steps = Int(round(total_time / dt))

set = se("system_ram.yaml")
set.segments = 2
set_values = [-50, -1.1, -1.1]

if !@isdefined s
    wing = RamAirWing(set)
    aero = BodyAerodynamics([wing])
    vsm_solver = Solver(aero; solver_type=NONLIN, atol=1e-8, rtol=1e-8)
    point_system = PointMassSystem(set, wing)
    s = RamAirKite(set, aero, vsm_solver, point_system)

    measure = Measurement()
end
s.set.abs_tol = 1e-5
s.set.rel_tol = 1e-3

if !ispath(joinpath(get_data_path(), "prob.bin"))
    KiteModels.init_sim!(s, measure)
end
measure.sphere_pos .= deg2rad.([50.0 50.0; 1.0 -1.0])
@time KiteModels.reinit!(s, measure)

logger = Logger(length(s.point_system.points), steps)

sys_state = KiteModels.SysState(s)
sys = s.sys
l = s.set.l_tether + 10
t = 0.
runtime = 0.
integ_runtime = 0.
try
    while t < total_time
        global t, runtime, integ_runtime
        PLOT && plot(s, t; zoom=false, front=false)
        global set_values = -s.set.drum_radius .* s.integrator[sys.winch_force]
        if (t < 3.0); set_values -= [0, 0, 5]; end
        steptime = @elapsed (t, integ_steptime) = next_step!(s, set_values; dt, vsm_interval)
        if (t > total_time/2); runtime += steptime; end
        if (t > total_time/2); integ_runtime += integ_steptime; end
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
p=plotx(logger.time_vec, 
        [logger.var_01_vec, logger.var_02_vec, logger.var_03_vec],
        [logger.var_04_vec, logger.var_05_vec],
        [logger.var_06_vec, logger.var_07_vec, logger.var_08_vec],
        [logger.var_09_vec, logger.var_10_vec, logger.var_11_vec, logger.var_12_vec],
        [logger.heading_vec]
        ;
    ylabels=["kite", "tether", "vsm", "twist", "heading"], 
    labels=[
        ["ω_b[1]", "ω_b[2]", "ω_b[3]"],
        ["vel[1]", "vel[2]"],
        ["force[3]", "kite moment[2]", "group moment[1]"],
        ["twist_angle[1]", "twist_angle[2]", "twist_angle[3]", "twist_angle[4]"],
        ["heading"]
        ],
    fig="Steering and heading MTK model")
display(p)

println("Times realtime: ", (total_time/2) / runtime)
println("Times realtime, just integrator: ", (total_time/2) / integ_runtime)

nothing