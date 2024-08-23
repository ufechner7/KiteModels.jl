using KiteModels, OrdinaryDiffEq, LinearAlgebra, Timers, SteadyStateDiffEq

using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using ControlPlots

update_settings()
set = se("system_3l.yaml")
set.abs_tol = 0.006
set.rel_tol = 0.01
steps = 110
dt = 1/set.sample_freq
tspan   = (0.0, dt)

logger = Logger(3*set.segments + 6, steps)

steering = [5,5,-30.0]

println("Running models")
mtk_kite::KPS4_3L = KPS4_3L(KCU(set))
pos = init_pos(mtk_kite)
simple_sys, _ = steady_state_model!(mtk_kite, pos; torque_control=false)
println("making steady state prob")
@time prob = SteadyStateProblem(ODEProblem(simple_sys, nothing, tspan))
println("solving steady state prob")
@time sol = solve(prob, DynamicSS(KenCarp4(autodiff=false)), abstol=0.1)
# @time sol = solve(prob, DynamicSS(Rodas5(autodiff=false)))
# println(sol)
pos0 = zeros(3, mtk_kite.num_A)
for i in 1:mtk_kite.num_A
    println(i)
    println("pos ", sol[simple_sys.pos_xz[2,i]])
    pos0[:,i] .= [sol[simple_sys.pos_xz[1,i]], pos[i][2], sol[simple_sys.pos_xz[2,i]]]
end
println("pos ", pos0)