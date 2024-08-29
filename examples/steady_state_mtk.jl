using KiteModels, OrdinaryDiffEqCore, OrdinaryDiffEqBDF, OrdinaryDiffEqSDIRK, LinearAlgebra, Timers
using Base: summarysize
tic()

using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using ControlPlots

set = deepcopy(load_settings("system_3l.yaml"))
set.abs_tol = 0.006
set.rel_tol = 0.01
steps = 110
dt = 1/set.sample_freq
tspan   = (0.0, dt)

logger = Logger(3*set.segments + 6, steps)

steering = [5,5,-30.0]

println("Running models")
s::KPS4_3L = KPS4_3L(KCU(set))
integrator = init_sim!(s; torque_control=true, mtk=true)

sys_state = KiteModels.SysState(s)
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
        steering = [5,5,-26.0] # left right middle
    end
    if i == 20
        steering = [0,10,-33]
    end
    if i == 40
        steering = [0,0,-10]
    end
    if i == 60
        steering = [0,0,-20]
    end

    if sys_state.heading > pi
        sys_state.heading -= 2*pi
    end
    sys_state.var_01 =  s.steering_pos[1]
    sys_state.var_02 =  s.steering_pos[2]
    sys_state.var_03 =  s.reel_out_speeds[1]
    sys_state.var_04 =  s.reel_out_speeds[2]

    total_new_time += @elapsed next_step!(s, integrator; set_values=steering)

    KiteModels.update_sys_state!(sys_state, s)
    if sys_state.heading > pi
        sys_state.heading -= 2*pi
    end
    # reltime = i*dt-dt
    # if mod(i, 5) == 1
    #     plot2d(s.pos, reltime; zoom=false, front=false, 
    #                             segments=set.segments, fig="side_view")            
    # end
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
