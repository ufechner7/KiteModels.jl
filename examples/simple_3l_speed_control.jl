using KiteModels, OrdinaryDiffEqCore, OrdinaryDiffEqBDF, OrdinaryDiffEqSDIRK, LinearAlgebra, Timers, Statistics
using Base: summarysize
tic()

using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using ControlPlots

set = deepcopy(load_settings("system_3l.yaml"))
# set.elevation = 71
dt = 0.05
total_time = 12.0

steps = Int(round(total_time / dt))
logger = Logger(3*set.segments + 6, steps)

if !@isdefined s; s = KPS4_3L(KCU(set)); end
s.set = update_settings()
s.set.abs_tol = 0.0006
s.set.rel_tol = 0.001
s.set.l_tether = 50.0
s.set.damping = 473
s.set.elevation = 85
println("init sim")
@time KiteModels.init_sim!(s; prn=true, torque_control=false)
# @time next_step!(s; set_values=[0.0, 0.0, 0.0], dt=2.0)
println("acc ", mean(norm.(s.integrator[s.simple_sys.force])))
sys_state = KiteModels.SysState(s)
if sys_state.heading > pi
    sys_state.heading -= 2*pi
end

println("stepping")
total_step_time = 0.0
toc()
for i in 1:steps
    time = (i-1) * dt
    @show time
    # println("acc ", norm(s.integrator[s.simple_sys.acc]))
    global total_step_time, sys_state, steering
    # steering = [0.0,0.0,1000.0] # left right middle
    if time < 3
        steering = [0.0,0.0,0.0] # left right middle
    elseif time < 4
        steering = [0,-0.9,-0.0]
    elseif time < 6
        steering = [-0.6,0.0,-0.0]
    elseif time < 20
        steering = [+0.6, +0.6, 0.0]
    end

    if sys_state.heading > pi
        sys_state.heading -= 2*pi
    end
    sys_state.var_01 =  rad2deg(s.flap_angle[1]) - 10
    sys_state.var_02 =  rad2deg(s.flap_angle[2]) - 10
    sys_state.var_03 =  s.reel_out_speeds[1]
    sys_state.var_04 =  s.reel_out_speeds[2]
    sys_state.var_05 =  s.reel_out_speeds[3]
    sys_state.var_06 =  rad2deg(s.integrator[s.simple_sys.seg_flap_angle[div(s.set.segments, 2)]] - s.integrator[s.simple_sys.aoa[div(s.set.segments, 2)]])
    sys_state.var_07 =  s.integrator[s.simple_sys.ram_force[div(s.set.segments, 2)]]
    sys_state.var_08 =  s.integrator[s.simple_sys.cl_seg[div(s.set.segments, 2)]]
    sys_state.var_09 =  s.integrator[s.simple_sys.cd_seg[div(s.set.segments, 2)]]
    # sys_state.var_09 =  norm(s.D_C + s.D_D)

    step_time = @elapsed next_step!(s; set_values=steering, dt=dt)
    if time > total_time/2
        total_step_time += step_time
    end

    KiteModels.update_sys_state!(sys_state, s)
    if sys_state.heading > pi
        sys_state.heading -= 2*pi
    end
    log!(logger, sys_state)
    # plot2d(s.pos, time; zoom=false, front=false, xlim=(-30, 30), ylim=(0, 60))
end

times_reltime = (total_time/2) / total_step_time
println("times realtime MTK model: ", times_reltime)
# println("avg steptime MTK model:   ", total_step_time/steps)

p=plotx(logger.time_vec, 
            [logger.var_01_vec,  logger.var_02_vec], 
            [logger.var_03_vec,  logger.var_04_vec, logger.var_05_vec], 
            rad2deg.(logger.heading_vec), 
            [logger.var_06_vec, logger.var_07_vec], 
            [logger.var_08_vec, logger.var_09_vec]; 
        ylabels=["Steering", "Reelout speed", "Heading [deg]", "Angle / Force", "Force"], 
        labels=[
            ["Steering Pos C", "Steering Pos D"], 
            ["v_ro left", "v_ro right", "v_ro middle"], 
            "Heading",
            ["Flap angle", "Ram Force"] ,
            ["Lift", "Drag"]], 
        fig="Steering and Heading MTK model")
display(p)