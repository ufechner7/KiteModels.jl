# plot the lift and drag coefficients as function of angle of attack

using Printf
using KiteModels, KitePodModels, KiteUtils

set = deepcopy(load_settings("system_v9.yaml"))

using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using ControlPlots

set.abs_tol=0.0006
set.rel_tol=0.00001

# the following values can be changed to match your interest
dt = 0.05
set.solver="DFBDF" # IDA or DFBDF
STEPS = 450
PLOT = false
PRINT = true
STATISTIC = false
DEPOWER = 0.525:-0.0025:0.3675
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
    integrator = KiteModels.init_sim!(kps4, stiffness_factor=0.1, prn=STATISTIC)
    set_depower_steering(kps4.kcu, rel_depower, 0.0)
    simulate(kps4, integrator, logger, steps)
end

CL = zeros(length(DEPOWER))
CD = zeros(length(DEPOWER))
AOA = zeros(length(DEPOWER))

elev = set.elevation
i = 1
for depower in DEPOWER
    global elev, i, kps4
    local cl, cd, aoa
    logger = Logger(set.segments + 5, STEPS)
    set.depower = 100*depower
    set.depower_gain = 10
    set.v_wind = 14
    kcu = KCU(set)
    kps4 = KPS4(kcu)
    cl, cd = sim_cl_cd(kps4, logger, depower)
    elev = rad2deg(logger.elevation_vec[end])
    set.elevation = elev
    aoa = kps4.alpha_2
    CL[i] = cl
    CD[i] = cd
    AOA[i] = aoa
    if PRINT
        println("Depower: $depower, CL $(round(cl, digits=2)), CD: $(round(cd, digits=2)), aoa: $(round(aoa, digits=2)), CL/CD: $(round(cl/cd, digits=2))")
        println("elevation: $(round((elev), digits=2))")
    end
    if depower in [DEPOWER[begin+1], DEPOWER[end]] && PLOT
        p = plot(logger.time_vec, rad2deg.(logger.elevation_vec), fig="depower: $depower")
        display(p)
    end
    i+=1
end

cl = zeros(length(AOA))
cd = zeros(length(AOA))
for (i, alpha) in pairs(AOA)
    global cl, cd
    cl[i] = kps4.calc_cl(alpha)
    cd[i] = kps4.calc_cd(alpha)
end

display(plot(AOA, [CL, cl], xlabel="AOA [deg]", ylabel="CL", labels=["CL","cl"], fig="CL vs AOA"))
display(plot(AOA, [CD, cd], xlabel="AOA [deg]", ylabel="CD", labels=["CD","cd"], fig="CD vs AOA"))
# AOA= 0:0.05:20
# calc_cd1 = KiteModels.Spline1D(se().alpha_cd, se().cd_list)
# plot(AOA, calc_cd1.(AOA), fig="calc_cd1", xlabel="AOA [deg]", ylabel="CD")
