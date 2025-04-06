using KiteModels

set = deepcopy(se())
# set.solver="IDA"
set.solver="DFBDF"

kps4::KPS4 = KPS4(KCU(set))
STEPS = 200

integrator = KiteModels.init_sim!(kps4; stiffness_factor=0.035, prn=false)
kps4.set.version = 2
kps4.stiffness_factor = 3
force = maximum(spring_forces(kps4; prn=false))

# next_step!(kps4, integrator, set_speed = 0, dt=0.05)
# function nsteps(kps4, integrator)
#     for i=1:STEPS
#         next_step!(kps4, integrator; set_speed = 0, dt=0.05)
#     end
# end
# nsteps(kps4, integrator)
# nothing
