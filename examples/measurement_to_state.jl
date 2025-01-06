using Revise, KiteModels, ModelingToolkit, ControlPlots, LinearAlgebra, Statistics

function angle_between_vectors(v1, v2)
    cos_theta = dot(v1, v2) / (norm(v1) * norm(v2))
    θ = acos(clamp(cos_theta, -1.0, 1.0))  # Clamp to handle numerical errors
    return rad2deg(θ)  # Angle in radians
end

dt = 0.001
total_time = 0.2
steps = Int(round(total_time / dt))

set = se("system_3l.yaml")
set.abs_tol = 1e-3
set.rel_tol = 1e-3
set.segments = 3
set.aero_surfaces = 2
logger = Logger(3*set.segments + 4, steps)
s = KPS4_3L(KCU(set))
s.measure.winch_torque = [-0.0, -0.0, 0.0]
s.measure.tether_acc = [0, 0, 0]
s.measure.tether_length = [52., 52., 49.]
s.measure.distance = 48.9
s.measure.elevation_left = deg2rad(80)
s.measure.elevation_right = deg2rad(80)
s.measure.azimuth_left = deg2rad(1)
s.measure.azimuth_right = deg2rad(-1)
# s.measure.distance_acc = s.measure.tether_acc[3]

prob, sol, ss, u0map = model!(s; real=true)
sys_state = KiteModels.SysState(s)

pos = sol[ss.pos]
pos = [[pos[j, i] for j in 1:3] for i in 1:s.i_C]

l = s.set.l_tether+10
plot2d(pos, 0.0; zoom=true, front=false, xlim=(-l/2, l/2), ylim=(0, l))

# next_step!(s)
@show sol.retcode
# for (u, e) in zip(sol.resid, equations(prob.f.sys))
#     println(u, "\t", e)
# end

t = 0
try
    while t < total_time
        global t
        t = next_step!(s; set_values=s.measure.winch_torque, dt)
        KiteModels.update_sys_state!(sys_state, s)
        sys_state.var_01 = s.integrator[ss.ω_p[1]]
        sys_state.var_02 = s.integrator[ss.ω_p[2]]
        sys_state.var_03 = s.integrator[ss.ω_p[3]]
        sys_state.var_04 = s.integrator[ss.kite_pos[1]]
        sys_state.var_05 = s.integrator[ss.kite_pos[2]]
        sys_state.var_06 = s.integrator[ss.kite_pos[3]]
        sys_state.var_07 = s.integrator[ss.trailing_edge_angle[1]]
        sys_state.var_08 = s.integrator[ss.trailing_edge_angle[2]]
        log!(logger, sys_state)
        local pos = [[sys_state.X[i], sys_state.Y[i], sys_state.Z[i]] for i in 1:s.i_C+1]
        plot2d(pos, t; zoom=false, front=false, xlim=(-l/2, l/2), ylim=(0, l))
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
            [logger.orient_vec],
            [logger.acc_vec],
            [logger.var_01_vec, logger.var_02_vec, logger.var_03_vec],
            [logger.var_04_vec, logger.var_05_vec, logger.var_06_vec],
            [logger.var_07_vec, logger.var_08_vec],
            ;
        ylabels=["orientation", "acc", "angular vel", "kite pos", "trailing edge angle"], 
        labels=[
            ["quaternion"],
            ["acc"],
            ["ω_p[1]", "ω_p[2]", "ω_p[3]"],
            ["x", "y", "z"],
            ["left", "right"],
            ],
        fig="Steering and heading MTK model")
display(p)

nothing