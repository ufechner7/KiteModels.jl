using KiteModels, OrdinaryDiffEq, LinearAlgebra, Timers, SteadyStateDiffEq, SymbolicIndexingInterface

using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using ControlPlots

update_settings()
set = se("system_3l.yaml")
set.abs_tol = 0.006
set.rel_tol = 0.01
steps = 50
dt = 1/set.sample_freq
tspan   = (0.0, dt)

logger = Logger(3*set.segments + 6, steps)

steering = [5,5,-30.0]

println("Running models")
s::KPS4_3L = KPS4_3L(KCU(set))
pos = init_pos(s)
simple_sys, _ = model!(s, pos; torque_control=true)
println("making steady state prob")
@time prob = SteadyStateProblem(ODEProblem(simple_sys, nothing, tspan))
println("solving steady state prob")
@time sol = solve(prob, DynamicSS(KenCarp4(autodiff=false); tspan=tspan), abstol=1e-6, reltol=1e-6)


# use optimized values
s.mtk = true
s.torque_control = true
# simple_sys, _ = model!(s, pos; torque_control=true)
s.prob = ODEProblem(simple_sys, nothing, tspan)
s.prob = remake(s.prob; u0=sol.u)
integrator = OrdinaryDiffEq.init(s.prob, KenCarp4(autodiff=false); dt, abstol=s.set.abs_tol, reltol=s.set.rel_tol, save_on=false)
s.set_values_idx = parameter_index(integrator.f, :set_values)
s.v_wind_gnd_idx = parameter_index(integrator.f, :v_wind_gnd)
s.v_wind_idx = parameter_index(integrator.f, :v_wind)
s.stiffness_factor_idx = parameter_index(integrator.f, :stiffness_factor)
s.get_pos = getu(integrator.sol, simple_sys.pos[:,:])
s.get_steering_pos = getu(integrator.sol, simple_sys.steering_pos)
s.get_line_acc = getu(integrator.sol, simple_sys.acc[:,s.num_E-2])
s.get_kite_vel = getu(integrator.sol, simple_sys.vel[:,s.num_A])
s.get_winch_forces = getu(integrator.sol, simple_sys.force[:,1:3])
s.get_L_C = getu(integrator.sol, simple_sys.L_C)
s.get_L_D = getu(integrator.sol, simple_sys.L_D)
s.get_D_C = getu(integrator.sol, simple_sys.D_C)
s.get_D_D = getu(integrator.sol, simple_sys.D_D)
s.get_tether_lengths = getu(integrator.sol, simple_sys.tether_length)
s.get_tether_speeds = getu(integrator.sol, simple_sys.tether_speed)
update_pos!(s, integrator)
# plot2d(s.pos, 0.0; zoom=false, front=false)




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
        steering = [-0.1,10,-33]
    end
    if i == 50
        steering = [-10.0,0.0,-20]
    end
    if i == 70
        steering = [0,0, -25]
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
