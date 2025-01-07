using Revise, KiteModels, ModelingToolkit, LinearAlgebra, Statistics

plot = true
if plot
    using ControlPlots
end

function angle_between_vectors(v1, v2)
    cos_theta = dot(v1, v2) / (norm(v1) * norm(v2))
    θ = acos(clamp(cos_theta, -1.0, 1.0))  # Clamp to handle numerical errors
    return rad2deg(θ)  # Angle in radians
end

dt = 0.01
total_time = 1.5
steps = Int(round(total_time / dt))

set = se("system_3l.yaml")
set.segments = 2
set.aero_surfaces = 6
logger = Logger(3*set.segments + 4, steps)
s = KPS4_3L(KCU(set))
s.measure.winch_torque = [-1, -1, -50.0]
s.measure.tether_acc = [0, 0, 0]
s.measure.tether_length = [51., 51., 49.]
s.measure.distance = 49.2
s.measure.elevation_left = deg2rad(70)
s.measure.elevation_right = deg2rad(70)
s.measure.azimuth_left = deg2rad(1)
s.measure.azimuth_right = deg2rad(-1)
# s.measure.distance_acc = s.measure.tether_acc[3]

prob, ss, u0map = model!(s)
sys_state = KiteModels.SysState(s)
l = s.set.l_tether + 10
t = 0
runtime = 0.0
try
    while t < total_time
        global t, runtime
        local pos = [[sys_state.X[i], sys_state.Y[i], sys_state.Z[i]] for i in 1:s.i_C+1]
        plot && plot2d(pos, t; zoom=false, front=false, xlim=(-l/2, l/2), ylim=(0, l), segments=10)
        set_values = -s.set.drum_radius * s.get_tether_forces(s.integrator)
        if t < 3; set_values[1] -= 10.0; end
        steptime = @elapsed t = next_step!(s; set_values, dt)
        if (t > total_time/2); runtime += steptime; end
        KiteModels.update_sys_state!(sys_state, s)
        sys_state.var_01 = s.integrator[ss.ω_b[1]]
        sys_state.var_02 = s.integrator[ss.ω_b[2]]
        sys_state.var_03 = s.integrator[ss.ω_b[3]]
        sys_state.var_04 = s.integrator[ss.tether_length[1]]
        sys_state.var_05 = s.integrator[ss.tether_length[2]]
        sys_state.var_06 = s.integrator[ss.tether_length[3]]
        sys_state.var_07 = s.integrator[ss.trailing_edge_angle[1]]
        sys_state.var_08 = s.integrator[ss.trailing_edge_angle[2]]
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
if plot 
    p=plotx(logger.time_vec, 
            [logger.orient_vec],
            [logger.acc_vec],
            [logger.var_01_vec, logger.var_02_vec, logger.var_03_vec],
            [logger.var_04_vec, logger.var_05_vec, logger.var_06_vec],
            [logger.var_07_vec, logger.var_08_vec],
            [logger.heading_vec],
            ;
        ylabels=["orientation", "acc", "angular vel", "tether length", "trailing edge angle", "heading"], 
        labels=[
            ["quaternion"],
            ["acc"],
            ["ω_b[1]", "ω_b[2]", "ω_b[3]"],
            ["left", "right", "middle"],
            ["left", "right"],
            ["heading_y"]
            ],
        fig="Steering and heading MTK model")
    display(p)
end

println("Times realtime: ", (total_time/2) / runtime)

nothing