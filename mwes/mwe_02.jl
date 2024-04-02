using KiteModels

kps4::KPS4 = KPS4(KCU(se()))
STEPS = 200

integrator = KiteModels.init_sim!(kps4; stiffness_factor=0.035, prn=false)
next_step!(kps4, integrator, v_ro = 0, dt=0.05)
@timev next_step!(kps4, integrator, v_ro = 0, dt=0.05)
function nsteps(kps4, integrator)
    for i=1:STEPS
        next_step!(kps4, integrator, v_ro = 0, dt=0.05)
    end
end
bytes = @allocated nsteps(kps4, integrator)
println("Bytes per step: $(Int64(round(bytes/STEPS)))")
nothing
# 3.125 ms, 2.452 MB with IDA
# 7.38  ms, 0.804 MB with DFBDF