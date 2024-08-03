using Revise, KiteModels, OrdinaryDiffEq, ControlPlots
# integrator = init_sim!(s; modeling_toolkit=true)

# s = KPS4_3L(KCU(se()))
# dt = 1/s.set.sample_freq
# print("using old method")

# integrator = KiteModels.init_sim!(s; stiffness_factor=0.5, prn=true, mtk=false)
# for i in 1:3
#     @time step!(integrator, dt, true)
# end

s = KPS4_3L(KCU(se()))
dt = 0.0001
integrator, simple_sys = KiteModels.init_sim!(s; stiffness_factor=0.3, prn=true, mtk=true)
# simple_sys, sys = model!(s, s.pos, s.vel)
for i in 1:10
    println("stepping...")
    @time next_step!(s, integrator)
    plot2d(s.pos, 10; zoom=false, front=false, segments=se().segments)
    for (p, v, a) in zip(s.pos, s.vel, s.acc)
        println("p $p v $v a $a ")
    end
    sleep(1)
end

"""
debugging:
integrator.sol(0.000000001; idxs=simple_sys.F_steering_c)
"""

nothing
