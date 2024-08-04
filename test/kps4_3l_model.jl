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
# # simple_sys, sys = model!(s, s.pos, s.vel)
# for i in 1:3
#     println("stepping...")
#     @time next_step!(s, integrator)
#     plot2d(s.pos, 10; zoom=false, front=false, segments=se().segments)
#     for (p, v, a) in zip(s.pos, s.vel, s.acc)
#         println("p $p v $v a $a ")
#     end
#     sleep(1)
# end

"""
debugging:
integrator.sol(0.000000001; idxs=simple_sys.F_steering_c)

lift forces 
s.forces[s.num_C] [90.77840952796961, 165.41850136730986, 342.11333043413924]
s.L_C [-8.497202522311054e-16, 155.80076139240956, 383.6297003775512]
s.dL_dα [-9.08473977288504e-16, -100.70528355285708, 608.742755393975]
s.dα 0.3416666666666666
"""

nothing
