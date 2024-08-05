using Revise, KiteModels, OrdinaryDiffEq, ControlPlots, SteadyStateDiffEq

s = KPS4_3L(KCU(se()))
dt = 1/s.set.sample_freq
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
    @elapsed update_pos!(s, integrator)
    plot2d(s.pos, i-1; zoom=false, front=false, segments=se().segments)
end

println("re-init")
@time integrator = OrdinaryDiffEq.init(prob, TRBDF2(autodiff=false); dt=dt, abstol=s.set.abs_tol, reltol=s.set.rel_tol)

println("stepping")
steps = 100
total_time = 0.0
for i in 1:steps
    global total_time += @elapsed OrdinaryDiffEq.step!(integrator, dt, true)
    global total_time += @elapsed update_pos!(s, integrator)
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

for (i,p) in enumerate(integrator.sol(integrator.sol.t[end]; idxs=simple_sys.pos))
    if s.pos[s.num_E-2][1] == p
        println(i, " ", p)
    end
end
"""

nothing
