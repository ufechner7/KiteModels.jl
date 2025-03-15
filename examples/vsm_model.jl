using Revise, KiteModels, ModelingToolkit, LinearAlgebra, Statistics, VortexStepMethod

PLOT = true
if PLOT
    using ControlPlots
end

dt = 1.0
total_time = 40.
steps = Int(round(total_time / dt))

set = se("system_3l.yaml")
set.segments = 6
set.aero_surfaces = 6
logger = Logger(3*set.segments + 4, steps)

# if !@isdefined(s); s = KPSQ(KCU(set)); end
wing = KiteWing("data/ram_air_kite_body.obj", "data/ram_air_kite_foil.dat"; mass=set.mass, crease_frac=0.9)
aero = BodyAerodynamics([wing])
solver = Solver()
s = KPSQ(set, wing, aero, solver)

s.measure.set_values = [-0.5, -0.5, -60.0]
s.measure.tether_length = [51., 51., 49.]
s.measure.tether_vel = [0.015, 0.015, 0.782]
s.measure.tether_acc = [0.18, 0.18, 4.12]
s.measure.sphere_pos[1, 1] = deg2rad(81.36)
s.measure.sphere_pos[1, 2] = deg2rad(81.36)
s.measure.sphere_pos[2, 1] = deg2rad(1)
s.measure.sphere_pos[2, 2] = deg2rad(-1)
s.measure.sphere_vel .= [0.13 0.13; 0 0]
s.measure.sphere_acc .= [0.09 0.09; 0 0]
s.set.abs_tol = 0.001
s.set.rel_tol = 0.001
s.set.segments = 2
# s.measure.distance_acc = s.measure.tether_acc[3]

sys, defaults, guesses = KiteModels.model!(s)
s.simple_sys = sys
@time s.prob = ODEProblem(sys, defaults, (0.0, 0.01); guesses)
solver = FBDF( # https://docs.sciml.ai/SciMLBenchmarksOutput/stable/#Results
    autodiff=ModelingToolkit.AutoFiniteDiff()
)
s.integrator = OrdinaryDiffEqCore.init(s.prob, solver; dt, abstol=s.set.abs_tol, reltol=s.set.rel_tol, save_on=false)
KiteModels.plot(s, 0.0)

# @time init_sim!(s; force_new_sys=false, force_new_pos=false, prn=true, ϵ=0.0, init=false)
# sys_state = KiteModels.SysState(s)
# sys = s.simple_sys
# l = s.set.l_tether + 10
# t = 0.
# runtime = 0.
# try
#     while t < total_time
#         global t, runtime
#         local pos = [[sys_state.X[i], sys_state.Y[i], sys_state.Z[i]] for i in 1:s.i_C+1]
#         PLOT && plot2d(pos, t; zoom=false, front=false, xlim=(-l/2, l/2), ylim=(0, l), segments=10)
#         # global set_values = -s.set.drum_radius * KiteModels.tether_force(s)
#         global set_values = s.measure.set_values
#         # if t < 1.0; set_values[2] -= 0.0; end
#         steptime = @elapsed t = next_step!(s; set_values, dt)
#         if (t > dt); runtime += steptime; end
#         KiteModels.update_sys_state!(sys_state, s)
#         sys_state.var_01 = s.get_α_b()[1]
#         sys_state.var_02 = s.get_α_b()[2]
#         sys_state.var_03 = s.get_α_b()[3]
#         sys_state.var_04 = s.get_wind_scale_gnd()
#         sys_state.var_05 = s.get_distance_acc()
#         sys_state.var_06 = s.get_distance()
#         sys_state.var_07 = s.get_trailing_edge_angle()[1]
#         sys_state.var_08 = s.get_trailing_edge_angle()[2]
#         sys_state.var_09 = s.get_tether_force()[1]
#         sys_state.var_10 = s.get_tether_force()[2]
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
# p=plotx(logger.time_vec, 
#         [logger.acc_vec],
#         [logger.var_01_vec, logger.var_02_vec, logger.var_03_vec],
#         [logger.var_04_vec],
#         [logger.var_05_vec],
#         [logger.var_07_vec, logger.var_08_vec],
#         [logger.var_09_vec, logger.var_10_vec],
#         [logger.heading_vec],
#         ;
#     ylabels=["acc", "α", "pos", "acc", "te angle", "force", "heading"], 
#     labels=[
#         ["acc"],
#         ["α_b[1]", "α_b[2]", "α_b[3]"],
#         ["wind speed"],
#         ["distance"],
#         ["left", "right"],
#         ["left winch", "right winch"],
#         ["heading_y"]
#         ],
#     fig="Steering and heading MTK model")
# display(p)

# println("Total runtime: ", runtime)
# println("Times realtime: ", (total_time) / runtime)

nothing