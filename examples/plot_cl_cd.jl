# plot the lift and drag coefficients as function of angle of attack

using Printf
using KiteModels, KitePodModels, KiteUtils

set = deepcopy(se())

set.abs_tol=0.0006
set.rel_tol=0.00001

# the following values can be changed to match your interest
dt = 0.05
set.solver="DFBDF" # IDA or DFBDF
STEPS = 450
PLOT = true
PRINT = false
STATISTIC = false
DEPOWER = 0.40:-0.01:0.21
# end of user parameter section #

if PLOT
    using Pkg
    if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
        using TestEnv; TestEnv.activate()
    end
    using ControlPlots
end

function simulate(kps4, integrator, logger, steps)
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
    simulate(kps4, integrator, logger, steps)
end

CL = zeros(length(DEPOWER))
CD = zeros(length(DEPOWER))
AOA = zeros(length(DEPOWER))

elev = set.elevation
i = 1
for depower in DEPOWER
    global elev, i
    logger = Logger(set.segments + 5, STEPS)
    set.depower = 100*depower
    set.depower_gain = 10
    set.v_wind = 12
    kcu = KCU(set)
    kps4 = KPS4(kcu)
    cl, cd = sim_cl_cd(kps4, logger, depower)
    elev = rad2deg(logger.elevation_vec[end])
    set.elevation = elev
    aoa = kps4.alpha_2
    CL[i] = cl
    CD[i] = cd
    AOA[i] = aoa
    println("Depower: $depower, CL $(round(cl, digits=2)), CD: $(round(cd, digits=2)), aoa: $(round(aoa, digits=2)), CL/CD: $(round(cl/cd, digits=2))")
    println("elevation: $(round((elev), digits=2))")
    if depower in [DEPOWER[begin+1], DEPOWER[end]] && PLOT
        p = plot(logger.time_vec, rad2deg.(logger.elevation_vec), fig="depower: $depower")
        display(p)
    end
    i+=1
end

display(plot(AOA, CL, xlabel="AOA [deg]", ylabel="CL", fig="CL vs AOA"))
display(plot(AOA, CD, xlabel="AOA [deg]", ylabel="CD", fig="CD vs AOA"))
nothing