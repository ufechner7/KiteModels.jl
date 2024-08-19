# plot the lift and drag as function of angle of attack

using Printf
using KiteModels, KitePodModels, KiteUtils

set = deepcopy(se())

set.abs_tol=0.0006
set.rel_tol=0.00001

# the following values can be changed to match your interest
dt = 0.05
set.solver="DFBDF" # IDA or DFBDF
STEPS = 400
PLOT = true
PRINT = false
STATISTIC = false
DEPOWER = 0.22:0.01:0.34
# end of user parameter section #

kcu::KCU = KCU(set)
kps4::KPS4 = KPS4(kcu)

if PLOT
    using Pkg
    if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
        using TestEnv; TestEnv.activate()
    end
    using ControlPlots
end

function simulate(integrator, logger, steps)
    iter = 0
    for i in 1:steps
        KiteModels.next_step!(kps4, integrator; set_speed=0, dt)
        sys_state = KiteModels.SysState(kps4)
        log!(logger, sys_state)
        iter += kps4.iter
    end
    KiteModels.cl_cd(kps4)
end

function sim_cl_cd(kps4::KPS4, logger, rel_depower; steps=STEPS)
    integrator = KiteModels.init_sim!(kps4, stiffness_factor=0.5, prn=STATISTIC)
    set_depower_steering(kps4.kcu, rel_depower, 0.0)
    simulate(integrator, logger, steps)
end

for depower in DEPOWER
    logger = Logger(set.segments + 5, STEPS)
    cl, cd = sim_cl_cd(kps4, logger, depower)
    println("Depower: $depower, CL, CD  [N]: $(round(cl, digits=2)), $(round(cd, digits=2))")
    if depower in [DEPOWER[begin], DEPOWER[end]] && PLOT
        p = plot(logger.time_vec, rad2deg.(logger.elevation_vec), fig="depower: $depower")
        display(p)
    end
end
