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
total_time = 10.0

steps = Int(round(total_time / dt))
logger = Logger(3*set.segments + 6, steps)

if !@isdefined s; s = KPS4_3L(KCU(set)); end
s.set = update_settings()
s.set.abs_tol = 0.0006
s.set.rel_tol = 0.001
s.set.l_tether = 50.0
# s.set.damping = 473
# s.set.elevation = 85
println("init sim")
@time KiteModels.init_sim!(s; prn=true, torque_control=false)
# @time next_step!(s; set_values=[0.0, 0.0, 0.0], dt=2.0)
println("vel ", mean(norm.(s.integrator[s.simple_sys.force])))
sys_state = KiteModels.SysState(s)

println("stepping")
total_step_time = 0.0
toc()
steering = [0.0, 0.0, 0.0]
start_tether = copy(s.tether_lengths[3])
amount = 0.6
sign = 1
for i in 1:steps
    time = (i-1) * dt
    @show time
    # println("vel ", norm(s.integrator[s.simple_sys.vel]))
    global total_step_time, sys_state, steering, sign
    # steering = [0.0,0.0,1000.0] # left right middle
    if s.tether_lengths[1] > s.tether_lengths[2] + 0.1
        sign = -1
    elseif s.tether_lengths[1] < s.tether_lengths[2] - 0.1
        sign = 1
    end
    steering[1] += sign * dt * amount
    steering[2] -= sign * dt * amount

    sys_state.var_01 =  rad2deg(s.flap_angle[1]) - 10
    sys_state.var_02 =  rad2deg(s.flap_angle[2]) - 10
    sys_state.var_03 =  s.tether_lengths[1]
    sys_state.var_04 =  s.tether_lengths[2]
    sys_state.var_05 =  s.tether_lengths[3]
    sys_state.var_06 =  rad2deg(s.integrator[s.simple_sys.seg_flap_angle[div(s.set.aero_surfaces, 2)]] - s.integrator[s.simple_sys.aoa[div(s.set.aero_surfaces, 2)]])
    sys_state.var_07 =  rad2deg(s.integrator[s.simple_sys.flap_vel[1]])
    sys_state.var_08 =  norm(s.D_C)
    sys_state.var_09 =  norm(s.D_D)
    sys_state.var_10 =  (s.integrator[s.simple_sys.vel[:, s.num_E-3]]) ⋅ s.e_z
    sys_state.var_11 =  norm(s.integrator[s.simple_sys.vel[:, s.num_E-3]] .- (s.integrator[s.simple_sys.vel[:, s.num_E-3]]) ⋅ s.e_z)
    # sys_state.var_09 =  norm(s.D_C + s.D_D)

    # @show s.integrator[s.simple_sys.aoa[div(s.set.aero_surfaces, 2)]]

    step_time = @elapsed next_step!(s; set_values=steering, dt=dt)
    if time > total_time/2
        total_step_time += step_time
    end

    KiteModels.update_sys_state!(sys_state, s)
    log!(logger, sys_state)
    l = s.set.l_tether+10
    # plot2d(s.pos, time; zoom=true, front=false, xlim=(-l/2, l/2), ylim=(0, l))
end

times_reltime = (total_time/2) / total_step_time
println("times realtime MTK model: ", times_reltime)
# println("avg steptime MTK model:   ", total_step_time/steps)

p=plotx(logger.time_vec, 
            [logger.var_01_vec,  logger.var_02_vec], 
            [logger.var_03_vec,  logger.var_04_vec], 
            rad2deg.(logger.heading_vec), 
            [logger.var_06_vec, logger.var_07_vec], 
            [logger.var_08_vec, logger.var_09_vec],
            [logger.var_10_vec, logger.var_11_vec]; 
        ylabels=["Steering", "Length", "heading [deg]", "Angle / Force", "Force", "Vel"], 
        labels=[
            ["Steering Pos C", "Steering Pos D"], 
            ["Left tether", "Right tether"], 
            "heading",
            ["Flap angle", "Flap vel"] ,
            ["Drag C", "Drag D"],
            ["Vel par", "Vel perp"]],
        fig="Steering and heading MTK model")
display(p)