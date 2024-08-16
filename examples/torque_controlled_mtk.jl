using KiteModels, OrdinaryDiffEq, LinearAlgebra, Timers
using Base: summarysize
tic()

using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using ControlPlots

update_settings()
set = se("system_3l.yaml")
set.abs_tol = 0.006
set.rel_tol = 0.01
steps = 150
dt = 1/set.sample_freq
tspan   = (0.0, dt)

logger = Logger(3*set.segments + 6, steps)

steering = [5,5,-30.0]

println("Running models")
if ! @isdefined mtk_kite; mtk_kite = KPS4_3L(KCU(set)); end
if ! @isdefined mtk_integrator
    mtk_integrator = KiteModels.init_sim!(mtk_kite; stiffness_factor=0.1, prn=false, mtk=true, torque_control=true)
else 
    mtk_integrator = KiteModels.reset_sim!(mtk_kite; stiffness_factor=1.0)
end

println("compiling")
total_new_time = 0.0
for i in 1:5
    global total_new_time += @elapsed next_step!(mtk_kite, mtk_integrator; set_values=steering)
end
sys_state = KiteModels.SysState(mtk_kite)
if sys_state.heading > pi
    sys_state.heading -= 2*pi
end
log!(logger, sys_state)

println("stepping")
total_old_time = 0.0
total_new_time = 0.0
toc()
for i in 1:steps
    global total_new_time, sys_state, steering
    if i == 1
        steering = [5,5,-30.0] # left right middle
    end
    if i == 20
        steering = [10,10,-30]
    end
    if i == 50
        steering = [0,10.0,-40]
    end

    if sys_state.heading > pi
        sys_state.heading -= 2*pi
    end
    sys_state.var_01 =  mtk_kite.steering_pos[1]
    sys_state.var_02 =  mtk_kite.steering_pos[2]
    sys_state.var_03 =  mtk_kite.reel_out_speeds[1]
    sys_state.var_04 =  mtk_kite.reel_out_speeds[2]

    total_new_time += @elapsed next_step!(mtk_kite, mtk_integrator; set_values=steering)

    sys_state = KiteModels.SysState(mtk_kite)
    if sys_state.heading > pi
        sys_state.heading -= 2*pi
    end
    log!(logger, sys_state)
end

new_time = (dt*steps) / total_new_time
println("times realtime MTK model: ", new_time)
println("avg steptime MTK model:   ", total_new_time/steps)

plotx(logger.time_vec, [logger.var_01_vec,  logger.var_02_vec], [logger.var_03_vec,  logger.var_04_vec], 
      rad2deg.(logger.heading_vec); 
      ylabels=["Steering", "Reelout speed", "Heading [deg]"], 
      labels=[["Steering Pos 1", "Steering Pos 2"], ["v_ro 1", "v_ro 2"], "Heading"], 
      fig="Steering and Heading MTK model")
