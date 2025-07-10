# SPDX-FileCopyrightText: 2025 Uwe Fechner
#
# SPDX-License-Identifier: MIT

using KiteModels
# test to check for bug in NLSolve package (or one of its dependencies)
# run it three times, if it always prints a force > 9000 newton no buggy package is in the dependencies

set = deepcopy(se())
set.solver="DFBDF"
set.v_wind = 50

kps4::KPS4 = KPS4(KCU(set))
STEPS = 200

integrator = KiteModels.init!(kps4; stiffness_factor=0.5, prn=false)

next_step!(kps4, integrator, set_speed = 0, dt=0.05)
function nsteps(kps4, integrator)
    for i=1:STEPS
        next_step!(kps4, integrator; set_speed = 0, dt=0.05)
    end
end
nsteps(kps4, integrator)

kps4.set.version = 2
kps4.stiffness_factor = 1
force = maximum(spring_forces(kps4; prn=false))
