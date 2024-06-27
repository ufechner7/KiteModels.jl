using Revise
using KiteModels
using ControlPlots
using LinearAlgebra, StaticArrays
update_settings()
kcu = KCU(se())
s = KPS4_3L(kcu)
integrator = KiteModels.init_sim!(s, stiffness_factor=0.04, prn=true, integrator_history=nothing)
for i in 1:200
# i = 1
    dt = 0.2
    plot2d(s.pos, i*dt; zoom=false, segments=s.set.segments)
    println(i)
    # println("connection lengths\t", s.l_connections)
    # println("force connection\t", s.forces[s.num_E-2])
    # println("speed \t", s.vel_connection)
    # for j in 4:s.num_A
    #     println("acc\t", (SVector(0, 0, -9.81) .+ s.forces[j] ./ s.masses[j]))
    # end
    @time KiteModels.next_step!(s, integrator, v_ro=[0.0,0.0,0.0], dt=dt)
    println("acc left\t", (SVector(0, 0, -9.81) .+ s.forces[s.num_E-2] ./ s.masses[s.num_E-2]))
    println("steering left\t", s.δ_left)
    println("acc right\t", (SVector(0, 0, -9.81) .+ s.forces[s.num_E-1] ./ s.masses[s.num_E-1]))
    println("steering right\t", s.δ_right)
    println("acc C\t", (SVector(0, 0, -9.81) .+ s.forces[s.num_C] ./ s.masses[s.num_C]))
    println("acc D\t", (SVector(0, 0, -9.81) .+ s.forces[s.num_D] ./ s.masses[s.num_D]))
    println("lift\t", s.lift_force)
end