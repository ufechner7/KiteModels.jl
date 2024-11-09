# plot the lift and drag coefficients as function of angle of attack

using Printf
using Pkg
using KiteModels, KitePodModels, KiteUtils, LinearAlgebra

set = deepcopy(load_settings("system.yaml"))

using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end

using ControlPlots
plt.close("all")

set.abs_tol=0.00006
set.rel_tol=0.000001
V_WIND = 14.5

# the following values can be changed to match your interest
dt = 0.05
set.solver="DFBDF" # IDA or DFBDF
STEPS = 1
PLOT = true
PRINT = true
STATISTIC = false
DEPOWER = 0.47:-0.005:0.355
# end of user parameter section #

elev = set.elevation
i = 1
set.v_wind = V_WIND # 25
logger::Logger = Logger(set.segments + 1, STEPS)

kcu::KCU = KCU(set)
kps3::KPS3 = KPS3(kcu)
integrator = KiteModels.init_sim!(kps3; delta=0.03, stiffness_factor=0.1, prn=STATISTIC)
lift, drag = lift_drag(kps3)
sys_state = KiteModels.SysState(kps3)
log!(logger, sys_state)
elev = rad2deg(logger.elevation_vec[end])
println("Lift: $lift, Drag: $drag, elev: $elev, Iterations: $(kps3.iter)")

q = QuatRotation(calc_orient_quat(kps3; viewer=false))
roll, pitch, yaw = rad2deg.(quat2euler(q))
println("--> orient_quat:       roll: ", roll, " pitch: ", pitch, "  yaw: ", yaw)
roll, pitch, yaw = rad2deg.(orient_euler(kps3))
println("--> orient_euler:      roll: ", roll, " pitch: ", pitch, " yaw:  ", yaw)
q = QuatRotation(calc_orient_quat(kps3; viewer=true))
roll, pitch, yaw = rad2deg.(quat2euler(q))
println("--> orient_quat (viewer): roll: ", roll, " pitch: ", pitch, "   yaw: ", yaw)

x, y, z = kite_ref_frame(kps3)
println("x:", x) # from trailing edge to leading edge in ENU reference frame
println("y:", y) # to the right looking in flight direction
println("z:", z) # down
azimuth = calc_azimuth(kps3)
println("azimuth: ", round(rad2deg(azimuth), digits = 2), "°")

# print alpha2, alpha3, alpha4
# println("alpha2, alpha3, alpha4: ", kps3.alpha_2, " ", kps3.alpha_3, " ", kps3.alpha_4)
println("heading: ", round(rad2deg(calc_heading(kps3)), digits = 2), "°")

ss = SysState(kps3)
println("AoA: ", rad2deg(ss.AoA), "°")
println("CL2: ", ss.CL2, " CD2: ", ss.CD2)
println("v_wind_gnd: ", ss.v_wind_gnd)
println("v_wind_kite: ", ss.v_wind_kite)

# output on main branch
# Lift: 1047.1795339611076, Drag: 281.39765463928745, elev: 72.77014, Iterations: 1552
# roll: -0.0 pitch: 5.755691707832273 yaw: 90.00000250447816
# x:[-0.9949585608372032, 0.0, 0.10028689952711267]
# y:[0.0, 1.0, 0.0]
# z:[-0.10028689952711267, -2.8981421002141166e-19, -0.9949585608372032]
# alpha2, alpha3, alpha4: 7.546125780476343 10.000000853773646 10.000000853773646

