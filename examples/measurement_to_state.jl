using Revise, KiteModels, ModelingToolkit, LinearAlgebra, Statistics

plot = true
if plot
    using ControlPlots
end

dt = 1.0
total_time = 50.0
steps = Int(round(total_time / dt))

set = se("system_3l.yaml")
set.segments = 6
set.aero_surfaces = 6
logger = Logger(3*set.segments + 4, steps)
if !@isdefined(s); s = KPSQ(KCU(set)); end
s.measure.set_values = [-0.5, -0.5, -60.0]
s.measure.tether_acc = [0, 0, 0]
s.measure.tether_length = [51., 51., 49.]
s.measure.sphere_pos[1, 1] = deg2rad(80)
s.measure.sphere_pos[1, 2] = deg2rad(80)
s.measure.sphere_pos[2, 1] = deg2rad(1)
s.measure.sphere_pos[2, 2] = deg2rad(-1)
s.measure.sphere_vel .= [0.2 0.2; 0 0]
s.set.abs_tol = 0.001
s.set.rel_tol = 0.0006
# s.measure.distance_acc = s.measure.tether_acc[3]

@time init_sim!(s; force_new_sys=true, prn=true, ϵ=0.0, init=false)
# @assert false
sys_state = KiteModels.SysState(s)
sys = s.simple_sys
l = s.set.l_tether + 10
t = 0
runtime = 0.0
try
    while t < total_time
        global t, runtime
        local pos = [[sys_state.X[i], sys_state.Y[i], sys_state.Z[i]] for i in 1:s.i_C+1]
        plot && plot2d(pos, t; zoom=false, front=false, xlim=(-l/2, l/2), ylim=(0, l), segments=10)
        # global set_values = -s.set.drum_radius * KiteModels.tether_force(s)
        global set_values = s.measure.set_values
        # if t < 1.0; set_values[2] -= 0.0; end
        steptime = @elapsed t = next_step!(s; set_values, dt)
        if (t > dt); runtime += steptime; end
        KiteModels.update_sys_state!(sys_state, s)
        sys_state.var_01 = s.get_α_b()[1]
        sys_state.var_02 = s.get_α_b()[2]
        sys_state.var_03 = s.get_α_b()[3]
        sys_state.var_04 = s.integrator[s.simple_sys.gust_factor]
        sys_state.var_05 = s.get_distance_acc()
        sys_state.var_06 = norm(s.get_acc()[:, 6])
        sys_state.var_07 = s.get_trailing_edge_angle()[1]
        sys_state.var_08 = s.get_trailing_edge_angle()[2]
        sys_state.var_09 = s.get_force()[:, s.i_A] ⋅ s.get_e_te_A()
        sys_state.var_10 = s.get_force()[:, s.i_B] ⋅ s.get_e_te_B()
        sys_state.var_11 = s.get_tether_force()[3]
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
            [logger.acc_vec],
            [logger.var_01_vec, logger.var_02_vec, logger.var_03_vec],
            [logger.var_04_vec],
            [logger.var_05_vec, logger.var_06_vec],
            [logger.var_07_vec, logger.var_08_vec],
            [logger.var_09_vec, logger.var_10_vec, logger.var_11_vec],
            [logger.heading_vec],
            ;
        ylabels=["acc", "α", "pos", "acc", "te angle", "force", "heading"], 
        labels=[
            ["acc"],
            ["α_b[1]", "α_b[2]", "α_b[3]"],
            ["gust factor"],
            ["distance", "middle tether"],
            ["left", "right"],
            ["left flap", "right flap", "middle tether"],
            ["heading_y"]
            ],
        fig="Steering and heading MTK model")
    display(p)
end

println("Total runtime: ", runtime)
println("Times realtime: ", (total_time) / runtime)

nothing