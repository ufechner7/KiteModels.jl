using KiteModels, LinearAlgebra

PLOT = true
if PLOT
    using ControlPlots
end

include("./plotting.jl")

dt = 0.05
total_time = 3.5
vsm_interval = 5
steps = Int(round(total_time / dt))

set = se("system_ram.yaml")
set.segments = 2
set_values = [-50, -1.1, -1.1]
set.quasi_static = true

# if !@isdefined s
    wing = RamAirWing(set)
    aero = BodyAerodynamics([wing])
    vsm_solver = Solver(aero; solver_type=NONLIN, atol=1e-8, rtol=1e-8)
    point_system = PointMassSystem(set, wing)
    s = RamAirKite(set, aero, vsm_solver, point_system)

    measure = Measurement()
# end
s.set.abs_tol = 1e-5
s.set.rel_tol = 1e-3

measure.sphere_pos .= deg2rad.([50.0 50.0; 1.0 -1.0])
if !ispath(joinpath(get_data_path(), "prob.bin"))
    KiteModels.init_sim!(s, measure)
end
@time KiteModels.reinit!(s, measure; reload=true)
sys = s.sys
s.integrator.ps[sys.steady] = true
next_step!(s; dt=10.0, vsm_interval=1)
s.integrator.ps[sys.steady] = false

logger = Logger(length(s.point_system.points), steps)

sys_state = KiteModels.SysState(s)
l = s.set.l_tether + 10
t = 0.
runtime = 0.
integ_runtime = 0.
try
    while t < total_time
        global t, runtime, integ_runtime
        PLOT && plot(s, t; zoom=true, front=false)
        global set_values = -s.set.drum_radius .* s.integrator[sys.winch_force]
        if (t < 2.0); set_values -= [0, 0, 5]; end
        steptime = @elapsed (t, integ_steptime) = next_step!(s, set_values; dt, vsm_interval)
        t -= 10.0
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

        sys_state.var_13 = clamp(s.integrator[sys.pulley_acc[1]], -100, 100)
        sys_state.var_14 = clamp(s.integrator[sys.pulley_acc[2]], -100, 100)

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
        [logger.var_13_vec, logger.var_14_vec],
        [logger.heading_vec]
        ;
    ylabels=["kite", "tether", "vsm", "twist", "pulley", "heading"], 
    labels=[
        ["ω_b[1]", "ω_b[2]", "ω_b[3]"],
        ["vel[1]", "vel[2]"],
        ["force[3]", "kite moment[2]", "group moment[1]"],
        ["twist_alpha[1]", "twist_alpha[2]", "twist_alpha[3]", "twist_alpha[4]"],
        ["pulley acc[1]", "pulley acc[2]"],
        ["heading"]
        ],
    fig="Steering and heading MTK model")
display(p)

println("Times realtime: ", (total_time/2) / runtime)
println("Times realtime, just integrator: ", (total_time/2) / integ_runtime)

@show norm(logger.var_13_vec[steps÷2:end])

# diffs = [norm(diff([u[i] for u in s.integrator.sol.u[1000:end]])) for i in eachindex(s.integrator.u)]
# names = unknowns(s.prob.f.sys)
# for (n, d) in zip(names, diffs)
#     println(round(d; digits=2), "\t", n)
# end


nothing