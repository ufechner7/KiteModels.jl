# plot the lift and drag coefficients as function of angle of attack

using Printf
using KiteModels, KitePodModels, KiteUtils

# angle_of_attack,          lift_coefficient,         drag_coefficient
# 2.000000000000000000e+00, 3.311587158641531303e-01, 1.110930319374049541e-01
# 3.000000000000000000e+00, 3.367252247881619698e-01, 9.189103773118127705e-02
# 4.000000000000000000e+00, 3.185656207459444111e-01, 9.895043587004535846e-02
# 5.000000000000000000e+00, 3.262272109321684987e-01, 9.438550669014250660e-02
# 6.000000000000000000e+00, 3.852419439968312598e-01, 1.001150540805576389e-01
# 7.000000000000000000e+00, 5.837620472869193833e-01, 1.387351148984416749e-01
# 8.000000000000000000e+00, 7.835204148162784321e-01, 1.867904340047445710e-01
# 9.000000000000000000e+00, 8.505108347128584878e-01, 2.011420208454205993e-01
# 1.000000000000000000e+01, 8.860003856047766746e-01, 2.052157022676944498e-01
# 1.100000000000000000e+01, 8.994221252043600456e-01, 2.015113631579626696e-01
# 1.200000000000000000e+01, 8.765915870833970169e-01, 1.873753922588730636e-01

set = deepcopy(load_settings("system_v9.yaml"))

alpha_cl=  [-180.0, -160.0, -90.0, -20.0, -10.0,  -5.0,  0.0, 2.0,     3.0,     4.0,     5.0,     6.0,     7.0,     8.0,     9.0,      10.0,   11.0,    12.0,       20.0,    40.0, 90.0, 160.0, 180.0]
cl_list=   [   0.0,    0.5,   0.0,  0.08, 0.125,  0.15,  0.2, 0.33115, 0.33672, 0.31856, 0.32623, 0.38524, 0.58376, 0.78352, 0.85051,  0.8860, 0.89942, 0.87659, 0.87659, 0.87659,  0.0,  -0.5,   0.0]
alpha_cd=  [-180.0, -170.0, -140.0, -90.0, -20.0, 0.0, 20.0, 90.0, 140.0, 170.0, 180.0]
cd_list=   [   0.5,    0.5,    0.5,   1.0,   0.2, 0.1,  0.2,  1.0,   0.5,   0.5,   0.5]

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
DEPOWER = 0.50:-0.005:0.34
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
AOA=-180:0.05:180
calc_cl1 = KiteModels.Spline1D(se().alpha_cl, se().cl_list)
plot(AOA, calc_cl1.(AOA), fig="calc_cl1", xlabel="AOA [deg]", ylabel="CL")
