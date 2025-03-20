using Revise, KiteModels, ModelingToolkit, LinearAlgebra, Statistics, VortexStepMethod
using OrdinaryDiffEqBDF, OrdinaryDiffEqCore, OrdinaryDiffEqRosenbrock
using ModelingToolkit: setp

PLOT = true
if PLOT
    using ControlPlots
end

dt = 0.005
total_time = 4.0
steps = Int(round(total_time / dt))

set = se("system_3l.yaml")
set.segments = 2
set_values = [-50, -2.1, -2.1]

"""
With aligned principal frame at t = 1.0
julia> s.integrator[sys.R_b_w]
3×3 Matrix{Float64}:
  0.957587  -0.00302636   0.288128
  0.003092   0.999995     0.000227269
 -0.288127   0.000673263  0.957592

With unaligned principal frame at t = 1.0
julia> s.integrator[sys.R_b_w]
3×3 Matrix{Float64}:
  0.948203   -0.0143791   0.31734
  0.0127404   0.999893    0.0072386
 -0.31741    -0.00282062  0.948284
"""

new_sys = false
if new_sys
    # if !@isdefined(s); s = KPSQ(KCU(set)); end
    wing = RamAirWing("data/ram_air_kite_body.obj", "data/ram_air_kite_foil.dat"; mass=set.mass, crease_frac=0.9)
    aero = BodyAerodynamics([wing])
    vsm_solver = Solver()
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
    @time s2.prob = ODEProblem(sys, defaults_, (0.0, 0.01); guesses=guesses_)
    s2.simple_sys = sys
    s = s2
else
    # s.wing = RamAirWing("data/ram_air_kite_body.obj", "data/ram_air_kite_foil.dat"; mass=set.mass, crease_frac=0.9)
    # s.aero = BodyAerodynamics([s.wing])
    # s.vsm_solver = Solver()

    VortexStepMethod.deform!(s.wing, zeros(s.wing.n_panels), zeros(s.wing.n_panels))
    VortexStepMethod.init!(s.aero)
    s.vsm_solver.sol.gamma_distribution .= 0.0
end
s.set.abs_tol = 1e-5
s.set.rel_tol = 1e-3

solver = FBDF()
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
        global set_values = -s.set.drum_radius .* s.integrator[sys.winch_force]
        # if t < 1.0; set_values[2] -= 0.0; end
        steptime = @elapsed t = next_step!(s; set_values, dt)
        if (t > dt); runtime += steptime; end
        KiteModels.update_sys_state!(sys_state, s)
        sys_state.var_01 = s.integrator[sys.α_b[3]]
        sys_state.var_02 = s.integrator.ps[sys.aero_kite_moment_b[3]]
        sys_state.var_03 = s.integrator[sys.moment_b[3]]

        sys_state.var_04 = s.integrator[sys.tether_vel[1]]
        sys_state.var_05 = s.integrator[sys.tether_vel[3]]

        sys_state.var_06 = norm(s.vsm_solver.sol.aero_force)
        sys_state.var_07 = s.vsm_solver.sol.aero_moments[2]
        sys_state.var_08 = sum(s.vsm_solver.sol.moment_distribution)

        sys_state.var_13 = s.integrator[sys.twist_angle[1]]
        sys_state.var_14 = s.integrator[sys.twist_angle[4]]

        sys_state.var_15 = norm(s.integrator[sys.acc[:, 9]])
        sys_state.var_16 = norm(s.integrator[sys.acc[:, 10]])
        println(
            "Va[1]: ", s.get_va_body(s.integrator)[1], 
            "\tTwist: ", sum(s.get_twist(s.integrator)), 
            "\tM[3]: ", s.vsm_solver.sol.aero_moments[3],
            "\tω_b: ", (s.integrator[sys.ω_b]),
            )
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
        [logger.var_13_vec, logger.var_14_vec],
        [logger.var_15_vec, logger.var_16_vec]
        ;
    ylabels=["z acc", "tether", "coefficients", "twist angle", "bridle"], 
    labels=[
        ["α_b[3]", "aero_kite_moment_b[3]", "moment_b[3]"],
        ["vel[1]", "vel[2]"],
        ["force", "torque[2]", "moment"],
        ["angle[1]", "angle[4]"],
        ["acc[9]", "acc[10]"]
        ],
    fig="Steering and heading MTK model")
display(p)

println("Total runtime: ", runtime)
println("Times realtime: ", (total_time) / runtime)

nothing