using KiteModels
using ControlPlots
using LinearAlgebra, StaticArrays
update_settings()
kcu = KCU(se())
s = KPS4_3L(kcu)
integrator = KiteModels.init_sim!(s, stiffness_factor=0.04, prn=true, integrator_history=nothing)

for i in 1:20
    # i = 1
    dt = 0.01
    plot2d(s.pos, i*dt; zoom=false, front=false, segments=s.set.segments)
    println(i)
    println("E\t", s.forces[s.num_E])
    println("C\t", s.forces[s.num_C])
    println("D\t", s.forces[s.num_D])
    println("Lift \t", s.lift_force)
    println("Drag \t", s.drag_force)
    # println("Pos \t", s.pos)
    @time KiteModels.next_step!(s, integrator, v_ro=[0.0,0.0,0.0], dt=dt)
end