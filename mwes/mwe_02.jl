# SPDX-FileCopyrightText: 2025 Uwe Fechner
#
# SPDX-License-Identifier: MIT

using KiteModels

set = deepcopy(se())
# set.solver="IDA"
set.solver="DFBDF"

kps4::KPS4 = KPS4(KCU(set))
STEPS = 200

integrator = KiteModels.init!(kps4; stiffness_factor=0.035, prn=false)
next_step!(kps4, integrator, set_speed = 0, dt=0.05)
@timev next_step!(kps4, integrator; set_speed = 0, dt=0.05)
function nsteps(kps4, integrator)
    for i=1:STEPS
        next_step!(kps4, integrator; set_speed = 0, dt=0.05)
    end
end
bytes = @allocated nsteps(kps4, integrator)
println("Bytes per step: $(Int64(round(bytes/STEPS)))")
nothing
# 1.76 ms, 1.101 MB with IDA   on desktop
# 4.70 ms, 0.793 MB with DFBDF on desktop
