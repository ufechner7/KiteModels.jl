include("../src/KiteModels.jl")
using .KiteModels
kps4_3l = KPS4_3L()
integrator = KiteModels.init_sim!(kps4_3l, stiffness_factor=0.04, prn=true, integrator_history=nothing)
