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
integrator = KiteModels.init_sim!(s; stiffness_factor=0.3, prn=true, mtk=true)
for i in 1:10
    @time step!(integrator, dt, true)
    plot2d(s.pos, 10; zoom=false, front=false, segments=se().segments)          
    # println(s.pos)
    sleep(1)
end

nothing
