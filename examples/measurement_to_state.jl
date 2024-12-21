using Revise, KiteModels, ModelingToolkit, ControlPlots, LinearAlgebra, Statistics

function angle_between_vectors(v1, v2)
    cos_theta = dot(v1, v2) / (norm(v1) * norm(v2))
    θ = acos(clamp(cos_theta, -1.0, 1.0))  # Clamp to handle numerical errors
    return rad2deg(θ)  # Angle in radians
end

dt = 0.01
total_time = 2.0
steps = Int(round(total_time / dt))

set = se("system_3l.yaml")
set.segments = 2
set.aero_surfaces = 2
logger = Logger(3*set.segments + 6, steps)
s = KPS4_3L(KCU(set))
s.measure.winch_torque = [-0.05, -0.05, -60]
s.measure.tether_length = [51, 51, 50]
s.measure.distance = 49
s.measure.elevation_left = deg2rad(80)
s.measure.elevation_right = deg2rad(80)
s.measure.azimuth_left = deg2rad(-1)
s.measure.azimuth_right = deg2rad(1)
s.measure.distance_acc = s.measure.tether_acc[3]

prob, sol, ss, u0map = model!(s; real=true)
sys_state = KiteModels.SysState(s)

pos = sol[ss.pos]
pos = [[pos[j, i] for j in 1:3] for i in 1:s.num_A]

l = s.set.l_tether+10
plot2d(pos, 0.0; zoom=true, front=false, xlim=(-l/2, l/2), ylim=(0, l))

# next_step!(s)
@show angle_between_vectors(pos[6], pos[6+3]-pos[6])
@show sol.retcode
# for (u, e) in zip(sol.resid, equations(prob.f.sys))
#     println(u, "\t", e)
# end

# t = 0
# try
#     while t < total_time
#         global t
#         plot2d(s.pos, t; zoom=true, front=false, xlim=(-l/2, l/2), ylim=(0, l))
#         t = next_step!(s; set_values=s.measure.winch_torque, dt)
#         KiteModels.update_sys_state!(sys_state, s)
#         sys_state.var_01 = s.get_tether_vels(s.integrator)[3]
#         sys_state.var_02 = mean(norm.([s.integrator[ss.acc[:, i]] for i in vcat(4:s.num_flap_C-1, s.num_flap_D+1:s.num_A)]))
#         sys_state.var_03 = angle_between_vectors(s.pos[6], s.pos[6+3]-s.pos[6])
#         sys_state.var_04 = norm(s.integrator[ss.vel])
#         log!(logger, sys_state)
#     end
# catch e
#     if isa(e, AssertionError)
#         @show t
#     else
#         rethrow(e)
#     end
# end
# p=plotx(logger.time_vec, 
#             [logger.var_01_vec, logger.var_04_vec],
#             [logger.var_02_vec],
#             [logger.var_03_vec]; 
#         ylabels=["vel","acc", "pos"], 
#         labels=[
#             ["tether_vel","vel"],
#             ["acc"],
#             ["angle"],],
#         fig="Steering and heading MTK model")
# display(p)
nothing