using KiteModels

kps4_3L::KPS4_3L = KPS4_3L(KCU(se()))

integrator = KiteModels.init_sim!(kps4_3L; stiffness_factor=0.035, prn=false)
nothing