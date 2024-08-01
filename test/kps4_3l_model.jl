using Revise, KiteModels
s = KPS4_3L(KCU(se()))
# integrator = init_sim!(s; modeling_toolkit=true)
KiteModels.model!(s)