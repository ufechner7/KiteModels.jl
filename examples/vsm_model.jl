using Revise, KiteModels, ModelingToolkit, LinearAlgebra, Statistics, VortexStepMethod
using OrdinaryDiffEqBDF, OrdinaryDiffEqCore

PLOT = true
if PLOT
    using ControlPlots
end

dt = 0.001
total_time = 0.01
steps = Int(round(total_time / dt))

set = se("system_3l.yaml")
set.segments = 2

new_sys = false
if new_sys
    # if !@isdefined(s); s = KPSQ(KCU(set)); end
    wing = RamAirWing("data/ram_air_kite_body.obj", "data/ram_air_kite_foil.dat"; mass=set.mass, crease_frac=0.9)
    aero = BodyAerodynamics([wing])
    solver = Solver()
    s = KPSQ(set, wing, aero, solver)
    s.measure.set_values = [-14.1, -1.3, -1.3]
    s.measure.tether_length = [51., 51., 49.]
    s.measure.tether_vel = [0.015, 0.015, 0.782]
    s.measure.tether_acc = [0.18, 0.18, 4.12]
    s.measure.sphere_pos[1, 1] = deg2rad(50.)
    s.measure.sphere_pos[1, 2] = deg2rad(50.)
    s.measure.sphere_pos[2, 1] = deg2rad(1)
    s.measure.sphere_pos[2, 2] = deg2rad(-1)
    s.measure.sphere_vel .= [0.13 0.13; 0 0]
    s.measure.sphere_acc .= [0.09 0.09; 0 0]
    s.set.damping = 473.0 * 2
    s.set.abs_tol = 1e-4
    s.set.rel_tol = 1e-2
    # s.measure.distance_acc = s.measure.tether_acc[3]

    sys, defaults_, guesses_ = KiteModels.model!(s)
    @time s.prob = ODEProblem(sys, defaults_, (0.0, 0.01); guesses=guesses_)
    s.simple_sys = sys
    solver = FBDF(autodiff=ModelingToolkit.AutoFiniteDiff())
end
s.vsm_solver.sol.gamma_distribution .= 0.0
s.integrator = OrdinaryDiffEqCore.init(s.prob, solver; dt, abstol=s.set.abs_tol, reltol=s.set.rel_tol, save_on=false)
KiteModels.generate_getters!(s)
logger = Logger(length(s.point_system.points), steps)

# @time init_sim!(s; force_new_sys=false, force_new_pos=false, prn=true, ϵ=0.0, init=false)
sys_state = KiteModels.SysState(s)
sys = s.simple_sys
l = s.set.l_tether + 10
t = 0.
runtime = 0.
try
    while t < total_time
        global t, runtime
        KiteModels.plot(s, t)
        global set_values = -s.set.drum_radius * s.integrator[sys.winch_force]
        # if t < 1.0; set_values[2] -= 0.0; end
        steptime = @elapsed t = next_step!(s; set_values, dt)
        if (t > dt); runtime += steptime; end
        KiteModels.update_sys_state!(sys_state, s)
        sys_state.var_01 = s.integrator[sys.kite_acc[1]]
        sys_state.var_02 = s.integrator[sys.kite_acc[2]]
        sys_state.var_03 = s.integrator[sys.kite_acc[3]]
        sys_state.var_04 = s.integrator[sys.tether_vel[1]]
        sys_state.var_05 = s.integrator[sys.tether_vel[3]]
        sys_state.var_06 = norm(s.vsm_solver.sol.force_coefficients)
        sys_state.var_07 = s.vsm_solver.sol.moment_coefficients[2]
        sys_state.var_08 = sum(s.vsm_solver.sol.moment_coefficient_distribution)
        sys_state.var_09 = s.integrator[sys.twist_α[1]]
        sys_state.var_10 = s.integrator[sys.twist_α[2]]
        sys_state.var_11 = s.integrator[sys.twist_α[3]]
        sys_state.var_12 = s.integrator[sys.twist_α[4]]
        sys_state.var_13 = s.integrator[sys.twist_angle[1]]
        sys_state.var_14 = s.integrator[sys.twist_angle[2]]
        sys_state.var_15 = s.integrator[sys.twist_angle[3]]
        sys_state.var_16 = s.integrator[sys.twist_angle[4]]
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
        [logger.var_13_vec, logger.var_14_vec, logger.var_15_vec, logger.var_16_vec],
        ;
    ylabels=["kite", "tether", "coefficients", "twist acc", "twist angle"], 
    labels=[
        ["acc[1]", "acc[2]", "acc[3]"],
        ["vel[1]", "vel[2]"],
        ["force", "torque[2]", "moment"],
        ["α[1]", "α[2]", "α[3]", "α[4]"],
        ["angle[1]", "angle[2]", "angle[3]", "angle[4]"],
        ],
    fig="Steering and heading MTK model")
display(p)

println("Total runtime: ", runtime)
println("Times realtime: ", (total_time) / runtime)

nothing