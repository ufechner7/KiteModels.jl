using Revise, KiteModels, ModelingToolkit, LinearAlgebra, Statistics, VortexStepMethod
using OrdinaryDiffEqBDF, OrdinaryDiffEqCore

PLOT = true
if PLOT
    using ControlPlots
end

dt = 0.01
total_time = 2.0
steps = Int(round(total_time / dt))

set = se("system_3l.yaml")
set.segments = 2

new_sys = false
if new_sys
    # if !@isdefined(s); s = KPSQ(KCU(set)); end
    wing = KiteWing("data/ram_air_kite_body.obj", "data/ram_air_kite_foil.dat"; mass=set.mass, crease_frac=0.9)
    aero = BodyAerodynamics([wing])
    solver = Solver()
    s = KPSQ(set, wing, aero, solver)
    s.measure.set_values = [-0.5, -0.5, -60.0]
    s.measure.tether_length = [51., 51., 49.]
    s.measure.tether_vel = [0.015, 0.015, 0.782]
    s.measure.tether_acc = [0.18, 0.18, 4.12]
    s.measure.sphere_pos[1, 1] = deg2rad(80.)
    s.measure.sphere_pos[1, 2] = deg2rad(80.)
    s.measure.sphere_pos[2, 1] = deg2rad(1)
    s.measure.sphere_pos[2, 2] = deg2rad(-1)
    s.measure.sphere_vel .= [0.13 0.13; 0 0]
    s.measure.sphere_acc .= [0.09 0.09; 0 0]
    s.set.abs_tol = 0.001
    s.set.rel_tol = 0.001
    # s.measure.distance_acc = s.measure.tether_acc[3]

    sys, defaults_, guesses_ = KiteModels.model!(s)
    s.simple_sys = sys
    @time s.prob = ODEProblem(sys, defaults_, (0.0, 0.01); guesses_)
    solver = FBDF( # https://docs.sciml.ai/SciMLBenchmarksOutput/stable/#Results
        autodiff=ModelingToolkit.AutoFiniteDiff()
    )
end
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
        sys_state.var_04 = s.integrator[sys.pulley_l0[1]]
        sys_state.var_05 = s.integrator[sys.pulley_l0[2]]
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
        [logger.var_04_vec, logger.var_05_vec]
        ;
    ylabels=["kite", "pulley"], 
    labels=[
        ["acc[1]", "acc[2]", "acc[3]"],
        ["l0[1]", "l0[2]"]
        ],
    fig="Steering and heading MTK model")
display(p)

println("Total runtime: ", runtime)
println("Times realtime: ", (total_time) / runtime)

nothing