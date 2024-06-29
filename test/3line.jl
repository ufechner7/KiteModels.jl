using KiteModels
using ControlPlots
using LinearAlgebra, StaticArrays
update_settings()
kcu = KCU(se())
s = KPS4_3L(kcu)
integrator = KiteModels.init_sim!(s, stiffness_factor=0.04, prn=true, integrator_history=nothing)

for i in 1:100
    # i = 1
    dt = 0.1
    plot2d(s.pos, i*dt; zoom=false, front=false, segments=s.set.segments)
    println(i)
    println("v_a ", s.v_a)
    println("v_kite ", s.v_kite)
    println("v_c ", s.vel[s.num_C])
    # println("Pos \t", s.pos)
    @time KiteModels.next_step!(s, integrator, v_ro=[0.0,0.0,0.0], dt=dt)
end