include("../src/KiteModels.jl")
using .KiteModels
kps4_3l = KPS4_3L()
integrator = KiteModels.init_sim!(kps4_3l, stiffness_factor=0.04, prn=true, integrator_history=nothing)
KiteModels.next_step!(kps4_3l, integrator, v_ro=[0.0,0.0,0.0], dt=0.2)