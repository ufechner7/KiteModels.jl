using KiteModels

kps4::KPS4 = KPS4(KCU(se()))

integrator = KiteModels.init_sim!(kps4; stiffness_factor=0.035, prn=false)
nothing