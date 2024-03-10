using KiteUtils, KiteModels, KitePodModels

SEGMENTS::Int64 = se().segments
kcu::KCU = KCU(se())
kps4::KPS4 = KPS4(kcu)

kps4.set.depower = 23.6
integrator = KiteModels.init_sim!(kps4; stiffness_factor=0.035, prn=false)
nothing
