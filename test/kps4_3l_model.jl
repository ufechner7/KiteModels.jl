using Revise, KiteModels, OrdinaryDiffEq
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
integrator = KiteModels.init_sim!(s; stiffness_factor=0.5, prn=true, mtk=true)
for i in 1:3
    @time step!(integrator, dt, true)
    println(integrator.u)
end

nothing

# There are 1420 highest order derivative variables and 1374 equations