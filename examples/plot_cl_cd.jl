# plot the lift and drag coefficients as function of angle of attack

using Printf
using KiteModels, KitePodModels, KiteUtils

set = deepcopy(load_settings("system_v9.yaml"))

using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using ControlPlots
plt.close("all")

set.abs_tol=0.00006
set.rel_tol=0.000001
V_WIND = 14

# the following values can be changed to match your interest
dt = 0.05
set.solver="DFBDF" # IDA or DFBDF
STEPS = 500
PLOT = true
PRINT = true
STATISTIC = false
DEPOWER = 0.45:-0.005:0.37
# DEPOWER = 0.41:-0.005:0.37
# end of user parameter section #

function set_tether_diameter!(se, d; c_spring_4mm = 614600, damping_4mm = 473)
    set.d_tether = d
    set.c_spring = c_spring_4mm * (d/4.0)^2
    set.damping = damping_4mm * (d/4.0)^2
end

set_tether_diameter!(set, set.d_tether)

if PLOT
    using Pkg
    if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
        using TestEnv; TestEnv.activate()
    end
    using ControlPlots
end

function simulate(kps4, integrator, logger, steps)
    iter = 0
    cl = 0.0
    cd = 0.0
    for i in 1:steps
        KiteModels.next_step!(kps4, integrator; set_speed=0, dt)
        sys_state = KiteModels.SysState(kps4)
        log!(logger, sys_state)
        iter += kps4.iter
        if i > steps - 50 # last 2.5s
            cl_, cd_ = KiteModels.cl_cd(kps4)
            cl += cl_
            cd += cd_
        end
    end
    return cl/50, cd/50
end


CL = zeros(length(DEPOWER)-2)
CD = zeros(length(DEPOWER)-2)
AOA = zeros(length(DEPOWER)-2)
DEP = zeros(length(DEPOWER)-2)

elev = set.elevation
i = 1
set.v_wind = V_WIND # 25
for depower in DEPOWER
    global elev, i, kps4
    local cl, cd, aoa, kcu
    depower == 0.41 && continue
    depower == 0.405 && continue
    logger = Logger(set.segments + 5, STEPS)
    DEP[i] = depower
    set.depower = 100*depower
    set.depower_gain = 5
    if i == 2
        set.v_wind = V_WIND
    end

    kcu = KCU(set)
    kps4 = KPS4(kcu)
    integrator = KiteModels.init_sim!(kps4; delta=0.03, stiffness_factor=0.05, prn=STATISTIC)
    if ! isnothing(integrator)
        set_depower_steering(kps4.kcu, depower, 0.0)
        cl, cd = simulate(kps4, integrator, logger, STEPS)
    else
        break
    end
    elev = rad2deg(logger.elevation_vec[end])
    if elev > 68
        set.v_wind = V_WIND - 2
    elseif elev > 64
        set.v_wind = V_WIND - 1
    end
    if elev > 50 && elev < 68
        set.elevation = elev
    end
    aoa = kps4.alpha_2
    CL[i] = cl
    CD[i] = cd
    AOA[i] = aoa
    if PRINT
        print("Depower: $depower, CL $(round(cl, digits=3)), CD: $(round(cd, digits=3)), aoa: $(round(aoa, digits=2)), CL/CD: $(round(cl/cd, digits=2))")
        println(", elevation: $(round((elev), digits=2))")
    end
    # if depower in [DEPOWER[begin+1], DEPOWER[end]] && PLOT
    if PLOT
        p = plot(logger.time_vec, rad2deg.(logger.elevation_vec), xlabel="time [s]", ylabel="elevation [°]", 
                 fig="depower: $depower")
        display(p)
        sleep(0.2)
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
display(plot(DEP, AOA, xlabel="Depower", ylabel="AOA [deg]", fig="AOA vs Depower"))
