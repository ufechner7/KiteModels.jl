using Revise, KiteModels, OrdinaryDiffEq, LinearAlgebra
using Base: summarysize

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
s = KPS4_3L(KCU(set))
integrator = KiteModels.init_sim!(s; stiffness_factor=0.1, prn=false, mtk=false, torque_control=true)
sys_state = KiteModels.SysState(s)
log!(logger, sys_state)

println("compiling")
total_new_time = 0.0
for i in 1:5
    global total_new_time += @elapsed next_step!(s, integrator; set_values=steering)
end
println("stepping")
total_old_time = 0.0
total_new_time = 0.0
steering_poss = [[],[]]
reel_out_speedss = [[],[]]
headings = []
for i in 1:steps
    global total_new_time, sys_state
    if i == 1
        global steering = [5,5,-30.0] # left right middle
    end
    if i == 20
        global steering = [10,10,-30]
    end
    if i == 50
        global steering = [0,10.0,-40]
    end

    if sys_state.heading > pi
        sys_state.heading -= 2*pi
    end
    sys_state.var_01 =  s.steering_pos[1]
    sys_state.var_02 =  s.steering_pos[2]
    sys_state.var_03 =  s.reel_out_speeds[1]
    sys_state.var_04 =  s.reel_out_speeds[2]
    total_new_time += @elapsed next_step!(s, integrator; set_values=steering)

    sys_state = KiteModels.SysState(s)
    if sys_state.heading > pi
        sys_state.heading -= 2*pi
    end
    log!(logger, sys_state)
end

new_time = (dt*steps) / total_new_time
println("times realtime new model: ", new_time)
println("avg steptime new model: ", total_new_time/steps)

plotx(logger.time_vec, [logger.var_01_vec,  logger.var_02_vec], [logger.var_03_vec,  logger.var_04_vec], rad2deg.(logger.heading_vec); 
      ylabels=["Steering", "Reelout speed", "Heading [deg]"], 
      labels=[["Steering Pos 1", "Steering Pos 2"], ["v_ro 1", "v_ro 2"], "Heading"], title="Steering and Heading")
