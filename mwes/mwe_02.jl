using KiteModels, Profile, StatProfilerHTML

kps4::KPS4 = KPS4(KCU(se()))

integrator = KiteModels.init_sim!(kps4; stiffness_factor=0.035, prn=false)
@timev next_step!(kps4, integrator, v_ro = 0, dt=0.05)
function nsteps(kps4, integrator)
    for i=1:10
        next_step!(kps4, integrator, v_ro = 0, dt=0.05)
    end
end
nsteps(kps4, integrator)
Profile.clear()
@profilehtml nsteps(kps4, integrator)
nothing