using Revise, KiteModels, OrdinaryDiffEq, ControlPlots, LinearAlgebra
using Base: summarysize
# using SteadyStateDiffEq

update_settings()
steps = 100
dt = 1/se().sample_freq
tspan   = (0.0, dt) 
plot2d([[0,0,0]], 0)

steering = [0,0,-0.5]


println("Running models")
if ! @isdefined s; s = KPS4_3L(KCU(se())); end;
if ! @isdefined integrator; integrator = KiteModels.init_sim!(s; stiffness_factor=0.1, prn=true, mtk=true, torque_control=true);
else integrator = KiteModels.reset_sim!(s); end;

println("compiling")
total_new_time = 0.0
for i in 1:2
    global total_new_time += @elapsed next_step!(s, integrator; set_values=steering)
end
println("stepping")
total_old_time = 0.0
total_new_time = 0.0
last_iter = 0
for i in 1:steps
    global last_iter
    if i==1
        global steering = [0,0,-1] # left right middle
    end
    if i==5
        global steering = [0,0,-1]
    end
    if i==15
        global steering = [0,1.0,0.0]
    end
    if i==20
        global steering = [1.0,0,0]
    end
    global total_new_time += @elapsed next_step!(s, integrator; set_values=steering)
    plot2d(s.pos, i-1; zoom=false, front=false, segments=se().segments)
    println(s.vel_kite ⋅ s.e_x)
    last_iter = integrator.iter
    readline()
    # println(s.vel_kite)
end
# plot2d(s1.pos, steps; zoom=false, front=false, segments=se().segments)

new_time = (dt*steps) / total_new_time
println("times realtime new model: ", new_time)
println("avg steptime new model: ", total_new_time/steps)
