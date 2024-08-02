using Revise, KiteModels
s = KPS4_3L(KCU(se()))
# integrator = init_sim!(s; modeling_toolkit=true)
simple_sys, sys = KiteModels.model!(s);
nothing

# There are 1420 highest order derivative variables and 1374 equations