using Revise, KiteModels, OrdinaryDiffEq, ControlPlots, SteadyStateDiffEq

s = KPS4_3L(KCU(se()))
dt = 0.2
tspan   = (0.0, dt) 
# integrator, simple_sys = KiteModels.init_sim!(s; stiffness_factor=0.1, prn=true, mtk=true)
y0, yd0 = KiteModels.find_steady_state!(s; stiffness_factor=0.1, prn=true, mtk=true)
simple_sys, sys = model!(s, y0, yd0)
@time prob = ODEProblem(simple_sys, nothing, tspan)
# steady_prob = SteadyStateDiffEq.SteadyStateProblem(prob)
# println("finding steady state")
# @time sol = solve(steady_prob, SSRootfind(), abstol=1e-6, reltol=1e-6)

# prob = ODEProblem(simple_sys, sol.u, tspan)

integrator = OrdinaryDiffEq.init(prob, TRBDF2(autodiff=false); dt=dt, abstol=s.set.abs_tol, reltol=s.set.rel_tol)

println("stepping")
steps = 10
total_time = 0.0
for i in 1:steps
    global total_time += @elapsed OrdinaryDiffEq.step!(integrator, dt, true)
    println(integrator.sol(integrator.sol.t[end]; idxs=simple_sys.acc[:,s.num_E]))
    @elapsed update_pos!(s, integrator)
    plot2d(s.pos, i-1; zoom=false, front=false, segments=se().segments)
end

println("re-init")
@time integrator = OrdinaryDiffEq.init(prob, TRBDF2(autodiff=false); dt=dt, abstol=s.set.abs_tol, reltol=s.set.rel_tol)

println("stepping")
steps = 10
total_time = 0.0
for i in 1:steps
    global total_time += @elapsed OrdinaryDiffEq.step!(integrator, dt, true)
    println(integrator.sol(integrator.sol.t[end]; idxs=simple_sys.acc[:,s.num_E]))
    @elapsed update_pos!(s, integrator)
    plot2d(s.pos, i-1; zoom=false, front=false, segments=se().segments)
end

println("total time ", total_time)
println("dt ", dt)
println("times realtime: ", (dt*steps) / total_time)



# # simple_sys, sys = model!(s, s.pos, s.vel)
# for i in 1:10
#     println("stepping...")
#     @time next_step!(s, integrator)
#     plot2d(s.pos, 10; zoom=false, front=false, segments=se().segments)
#     for (p, v, a) in zip(s.pos, s.vel, s.acc)
#         println("p $p v $v a $a ")
#     end
#     # sleep(1)
# end

"""
debugging:
integrator.sol(0.000000001; idxs=simple_sys.F_steering_c)

lift forces are good

spring forces 10
spring_force [-0.0, -90.81354570975981, -0.0]
half_drag_force [0.16832563706628656, 0.0, 0.0]
spring forces 11
spring_force [1.2037143896827731, -1.2512766659834529, -0.38775322658779593]
half_drag_force [0.047756472684599656, 0.04191603077077581, 0.012989194651326289]
spring forces 8
spring_force [99.75600270965079, 72.9483244909897, 337.30094388204435]
half_drag_force [0.3654202337395544, -0.02232848783231474, -0.10324322147013475]
"""

nothing
