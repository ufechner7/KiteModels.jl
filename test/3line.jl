using Revise
using KiteModels
using ControlPlots
using LinearAlgebra
update_settings()
kcu = KCU(se())
s = KPS4_3L(kcu)
integrator = KiteModels.init_sim!(s, stiffness_factor=0.04, prn=true, integrator_history=nothing)
# for i in 1:200
i = 1
    dt = 0.2
    plot2d(s.pos, i*dt; zoom=false, segments=s.set.segments)
    println(i)
    println("line_length\t", s.l_tethers)

    @profview KiteModels.next_step!(s, integrator, v_ro=[0.0,0.0,0.0], dt=dt)
# end