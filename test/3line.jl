using KiteModels
update_settings()
kcu = KCU(se())
kps4_3l = KPS4_3L(kcu)
integrator = KiteModels.init_sim!(kps4_3l, stiffness_factor=0.04, prn=true, integrator_history=nothing)
println("getting next step")
KiteModels.next_step!(kps4_3l, integrator, v_ro=[0.0,0.0,0.0], dt=0.2)