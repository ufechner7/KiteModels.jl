using Revise, KiteModels, OrdinaryDiffEq, ControlPlots, LinearAlgebra
using Base: summarysize
# using SteadyStateDiffEq

update_settings()
steps = 50
dt = 1/se().sample_freq
tspan   = (0.0, dt) 
plot2d([[0,0,0]], 0)

steering = [0,0,-0.5]


println("Running models")
s2 = KPS4_3L(KCU(se()))
s1 = KPS4_3L(KCU(se()))
integrator2, simple_sys = KiteModels.init_sim!(s2; stiffness_factor=0.1, prn=true, mtk=true)
integrator1 = KiteModels.init_sim!(s1; stiffness_factor=0.1, prn=true, mtk=false)
println("compiling")
total_old_time = 0.0
total_new_time = 0.0
for i in 1:2
    global total_new_time += @elapsed next_step!(s2, integrator2; set_values=steering)
    global total_old_time += @elapsed next_step!(s1, integrator1; set_values=steering)
end
println("stepping")
total_old_time = 0.0
total_new_time = 0.0
for i in 1:steps
    if i==1
        global steering = [0,0,0.3] # left right middle
    end
    if i==5
        global steering = [0,0,-0.3]
    end
    if i==15
        global steering = [0,0.3,0.0]
    end
    if i==20
        global steering = [0,0,0]
    end
    println(s2.vel_kite)
    global total_new_time += @elapsed next_step!(s2, integrator2; set_values=steering)
    global total_old_time += @elapsed next_step!(s1, integrator1; set_values=steering)
    plot2d(s2.pos, i-1; zoom=false, front=false, segments=se().segments)
    # plot2d(vcat(s1.pos,s2.pos), i-1; zoom=false, front=false, segments=se().segments)
    # println(norm(s1.pos-s2.pos)/3.0)
    # for (i,u) in enumerate(integrator2.u)
    #     println("i ", i, ", ", u)
    # end
    # println(s2.vel_kite)
    println(s1.vel_kite)
    # println(integrator2.sol[simple_sys.lengths[:]][end])
    # println(integrator2.sol[simple_sys.reel_out_speed[:]][end])
end
# plot2d(s1.pos, steps; zoom=false, front=false, segments=se().segments)

old_time = (dt*steps) / total_old_time
new_time = (dt*steps) / total_new_time
println("times realtime old model: ", old_time)
println("times realtime new model: ", new_time)
println("avg steptime new model: ", total_new_time/steps)
# println("old pos ", old_pos)
# println("new pos ", new_pos)
println("times faster new model: ", total_old_time/total_new_time)

"""

times realtime old model: 3.429651068617265
times realtime new model: 98.6601089002389
times faster new model: 28.766806571962952

debugging:

for (i,p) in enumerate(integrator.sol(integrator.sol.t[end]; idxs=simple_sys.pos))
    if s.pos[s.num_E-2][1] == p
        println(i, " ", p)
    end
end

steady_prob = SteadyStateDiffEq.SteadyStateProblem(prob)
@time sol = solve(steady_prob, SSRootfind(), abstol=1e-6, reltol=1e-6)
prob = ODEProblem(simple_sys, sol.u, tspan)
"""

nothing
