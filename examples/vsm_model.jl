using Revise, KiteModels, ModelingToolkit, LinearAlgebra, Statistics, VortexStepMethod
using OrdinaryDiffEqBDF, OrdinaryDiffEqCore, OrdinaryDiffEqRosenbrock, Serialization
using ModelingToolkit: setp

PLOT = true
if PLOT
    using ControlPlots
end

dt = 0.5
total_time = 6.5
steps = Int(round(total_time / dt))

set = se("system_3l.yaml")
set.segments = 2
set_values = [-50, -1.1, -1.1]

new_sys = 3
if new_sys == 1
    # if !@isdefined(s); s = KPSQ(KCU(set)); end
    wing = RamAirWing("data/ram_air_kite_body.obj", "data/ram_air_kite_foil.dat"; mass=set.mass, crease_frac=0.82, align_to_principal=true)
    aero = BodyAerodynamics([wing])
    vsm_solver = Solver(aero; solver_type=NONLIN, atol=1e-8, rtol=1e-8)
    s2 = KPSQ(set, wing, aero, vsm_solver)
    s2.measure.set_values = set_values
    s2.measure.tether_length = [51., 51., 49.]
    s2.measure.tether_vel = [0.015, 0.015, 0.782]
    s2.measure.tether_acc = [0.18, 0.18, 4.12]
    s2.measure.sphere_pos[1, 1] = deg2rad(80.)
    s2.measure.sphere_pos[1, 2] = deg2rad(80.)
    s2.measure.sphere_pos[2, 1] = deg2rad(1)
    s2.measure.sphere_pos[2, 2] = deg2rad(-1)
    s2.measure.sphere_vel .= [0.13 0.13; 0 0]
    s2.measure.sphere_acc .= [0.09 0.09; 0 0]
    s2.set.damping = 473.0 * 2
    # s.measure.distance_acc = s.measure.tether_acc[3]

    sys, defaults_, guesses_ = KiteModels.model!(s2)
    @info "Creating problem"
    @time s2.prob = ODEProblem(sys, defaults_, (0.0, 0.01); guesses=guesses_)
    serialize("data/kite.bin", s2.prob)
    # @info "Deserializing problem"
    # @time s2.prob = deserialize("data/kite.bin")
    s2.simple_sys = s2.prob.f.sys
    s = s2
elseif new_sys == 2
    wing = RamAirWing("data/ram_air_kite_body.obj", "data/ram_air_kite_foil.dat"; mass=set.mass, crease_frac=0.82, align_to_principal=true)
    aero = BodyAerodynamics([wing])
    vsm_solver = Solver(aero)
    s = KPSQ(set, wing, aero, vsm_solver)
    s.prob = deserialize("data/kite.bin")
    s.simple_sys = s.prob.f.sys
    s.point_system = KiteModels.PointMassSystem(s, s.wing)
elseif new_sys == 3
    # VortexStepMethod.deform!(s.wing, zeros(s.wing.n_panels), zeros(s.wing.n_panels))
    # VortexStepMethod.init!(s.aero)
    # s.vsm_solver.sol.gamma_distribution .= 0.0
end
s.set.abs_tol = 1e-5
s.set.rel_tol = 1e-3

solver = FBDF()
y = s.get_y(s.prob)
jac, x = VortexStepMethod.linearize(
    s.vsm_solver, 
    s.aero, 
    y;
    theta_idxs=1:4,
    va_idxs=5:7,
    omega_idxs=8:10,
    moment_frac=s.bridle_fracs[s.point_system.groups[1].fixed_index])
s.set_vsm(s.prob, [x, y, jac])
@info "Initializing integrator"
@time s.integrator = OrdinaryDiffEqCore.init(s.prob, solver; dt, abstol=s.set.abs_tol, reltol=s.set.rel_tol, save_on=false)
KiteModels.generate_getters!(s)
logger = Logger(length(s.point_system.points), steps)

# @time init_sim!(s; force_new_sys=false, force_new_pos=false, prn=true, ϵ=0.0, init=false)
sys_state = KiteModels.SysState(s)
sys = s.simple_sys
l = s.set.l_tether + 10
t = 0.
runtime = 0.
integ_runtime = 0.
try
    while t < total_time
        global t, runtime, integ_runtime
        KiteModels.plot(s, t; zoom=false, front=true)
        global set_values = -s.set.drum_radius .* s.integrator[sys.winch_force] - [0, 0, 5]
        # if t < 1.0; set_values[2] -= 0.0; end
        vsm_interval = 1
        steptime = @elapsed (t, integ_steptime) = next_step!(s; set_values, dt, vsm_interval)
        if (t > total_time/2); runtime += steptime; end
        if (t > total_time/2); integ_runtime += integ_steptime; end
        KiteModels.update_sys_state!(sys_state, s)
        sys_state.var_01 = s.integrator[sys.ω_b[1]]
        # sys_state.var_02 = s.integrator[sys.ω_b[2]]
        sys_state.var_03 = s.integrator[sys.ω_b[3]]

        sys_state.var_04 = s.integrator[sys.tether_vel[1]]
        sys_state.var_05 = s.integrator[sys.tether_vel[3]]

        sys_state.var_06 = s.integrator[sys.aero_force_b[3]]
        sys_state.var_07 = s.integrator[sys.aero_moment_b[2]]
        sys_state.var_08 = s.integrator[sys.group_aero_moment[1]]

        sys_state.var_09 = 0.01s.integrator[sys.group_aero_moment[4]]
        sys_state.var_10 = 0.01s.integrator[sys.group_tether_moment[4]]
        sys_state.var_11 = s.integrator[sys.twist_α[4]]
        sys_state.var_12 = s.integrator[sys.twist_angle[4]]

        sys_state.var_13 = norm(s.integrator[sys.spring_force[3]])
        sys_state.var_14 = norm(s.integrator[sys.spring_force[7]])
        # println(
        #     "\tPos: ", s.integrator[sys.pos[:, 17]],
        #     )
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
        [logger.var_13_vec, logger.var_14_vec]
        ;
    ylabels=["z acc", "tether", "vsm", "twist", "bridle"], 
    labels=[
        ["ω_b[1]", "ω_b[2]", "ω_b[3]"],
        ["vel[1]", "vel[2]"],
        ["force[3]", "kite moment[2]", "group moment[1]"],
        ["aero_moment[4]", "tether_moment[4]", "twist_α[4]", "twist_angle[4]"],
        ["acc[9]", "acc[10]"]
        ],
    fig="Steering and heading MTK model")
display(p)

println("Times realtime: ", (total_time/2) / runtime)
println("Times realtime, just integrator: ", (total_time/2) / integ_runtime)

nothing