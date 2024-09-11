using KiteModels, OrdinaryDiffEqCore, OrdinaryDiffEqBDF, OrdinaryDiffEqSDIRK, LinearAlgebra, Timers
using Base: summarysize
tic()

using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using ControlPlots

set = deepcopy(load_settings("system_3l.yaml"))
# set.elevation = 71
dt = 0.001
total_time = 1.0

steps = Int(round(total_time / dt))
logger = Logger(3*set.segments + 6, steps)

if !@isdefined s; s = KPS4_3L(KCU(set)); end
s.set = update_settings()
s.set.abs_tol = 0.006
s.set.rel_tol = 0.01
s.set.l_tether = 50.1
s.set.damping *= 1
println("init sim")
integrator = KiteModels.init_sim!(s; prn=true, torque_control=true, stiffness_factor=1.0)
println("acc ", norm(integrator[s.simple_sys.acc]))
sys_state = KiteModels.SysState(s)
if sys_state.heading > pi
    sys_state.heading -= 2*pi
end
log!(logger, sys_state)

println("stepping")
total_step_time = 0.0
toc()
for i in 1:steps
    time = (i-1) * dt
    @show time
    # println("acc ", norm(integrator[s.simple_sys.acc]))
    global total_step_time, sys_state, steering
    if time < 0.5
        steering = [20,10,-200.0] # left right middle
    elseif time < 1.0
        steering = [20,10,-100]
    end
    # if i == 40
    #     steering = [0,0,-20]
    # end
    # if i == 60
    #     steering = [0,0,-30]
    # end

    if sys_state.heading > pi
        sys_state.heading -= 2*pi
    end
    sys_state.var_01 =  rad2deg(s.flap_angle[1])
    sys_state.var_02 =  rad2deg(s.flap_angle[2])
    sys_state.var_03 =  s.reel_out_speeds[1]
    sys_state.var_04 =  s.reel_out_speeds[2]
    sys_state.var_05 =  s.reel_out_speeds[3]
    sys_state.var_06 =  norm((integrator[s.simple_sys.acc[:, 6]] ⋅ normalize(s.pos[6])) * normalize(s.pos[6]))
    sys_state.var_07 =  norm(integrator[s.simple_sys.acc[:, 9]] .- (integrator[s.simple_sys.acc[:, 6]] ⋅ normalize(s.pos[6])) * normalize(s.pos[6]))

    step_time = @elapsed next_step!(s, integrator; set_values=steering, dt=dt)
    if time > 0.5
        total_step_time += step_time
    end

    KiteModels.update_sys_state!(sys_state, s)
    if sys_state.heading > pi
        sys_state.heading -= 2*pi
    end
    log!(logger, sys_state)
    # plot2d(s.pos, time; zoom=false, front=false, xlim=(0, 100), ylim=(0, 100))
end

times_reltime = (total_time - 0.5) / total_step_time
println("times realtime MTK model: ", times_reltime)
# println("avg steptime MTK model:   ", total_step_time/steps)

p=plotx(logger.time_vec, [logger.var_01_vec,  logger.var_02_vec], [logger.var_03_vec,  logger.var_04_vec, logger.var_05_vec], 
        rad2deg.(logger.heading_vec), [logger.var_06_vec, logger.var_07_vec]; 
        ylabels=["Steering", "Reelout speed", "Heading [deg]", "Acc"], 
        labels=[["Steering Pos C", "Steering Pos D"], ["v_ro left", "v_ro right", "v_ro middle"], "Heading", ["middle tether", "perp middle tether"]], 
        fig="Steering and Heading MTK model")
display(p)
