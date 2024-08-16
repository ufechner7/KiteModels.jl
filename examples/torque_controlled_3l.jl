using Revise, KiteModels, OrdinaryDiffEq, LinearAlgebra, Plots
using Base: summarysize
# using SteadyStateDiffEq

update_settings()
set = se("system_3l.yaml")
set.abs_tol = 0.006
set.rel_tol = 0.01
steps = 150
dt = 1/set.sample_freq
tspan   = (0.0, dt) 
# plot2d([[0,0,0]], 0)

steering = [5,5,-30.0]


println("Running models")
s = KPS4_3L(KCU(set))
integrator = KiteModels.init_sim!(s; stiffness_factor=0.1, prn=false, mtk=false, torque_control=true)

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
    if i == 1
        global steering = [5,5,-30.0] # left right middle
    end
    if i == 20
        global steering = [10,10,-30]
    end
    if i == 50
        global steering = [0,10.0,-40]
    end
    # if i == 10
    #     global steering = [-0.5,10.5,-0]
    # end
    # if i == 20
    #     global steering = [10.5,10.5,13]
    # end
    # println(s.steering_pos, "\t tether_lengths \t", s.reel_out_speeds)
    heading = calc_heading(s)
    # println(s.steering_pos)
    if heading > pi
        heading -= 2*pi
    end
    push!(headings, heading)
    push!(steering_poss[1], s.steering_pos[1])
    push!(steering_poss[2], s.steering_pos[2])
    push!(reel_out_speedss[1], s.reel_out_speeds[1])
    push!(reel_out_speedss[2], s.reel_out_speeds[2])
    # println(norm.(s.winch_forces))
    global total_new_time += @elapsed next_step!(s, integrator; set_values=steering)
    # plot2d(s.pos, i-1; zoom=false, front=true, segments=set.segments)
end
# plot2d(s1.pos, steps; zoom=false, front=false, segments=set.segments)

new_time = (dt*steps) / total_new_time
println("times realtime new model: ", new_time)
println("avg steptime new model: ", total_new_time/steps)

Plots.plot(1:steps, [steering_poss, reel_out_speedss, headings], 
    label=["Steering Pos 1" "Steering Pos 2" "Reel Out Speed 1" "Reel Out Speed 2" "Heading"],
    linewidth=3)

# 0.02147723333217222 rad with 500 damping
# damping was too much...

# 0.14807810376419894 with 50 damping
