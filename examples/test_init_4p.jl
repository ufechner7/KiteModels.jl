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
set.v_reel_out = 1.0 # initial reel-out speed [m/s]
STEPS = 1
PLOT = true
PRINT = true
STATISTIC = false
DEPOWER = 0.47:-0.005:0.355
UPWIND_DIR = -pi/2 +deg2rad(10)
# end of user parameter section #

elev = set.elevation
i = 1
set.v_wind = V_WIND # 25
logger::Logger = Logger(set.segments + 5, STEPS)

kcu::KCU = KCU(set)
kps4::KPS4 = KPS4(kcu)
integrator = KiteModels.init_sim!(kps4; delta=0.03, stiffness_factor=0.01, upwind_dir=UPWIND_DIR, prn=STATISTIC)
for i=1:80
    KiteModels.next_step!(kps4, integrator; set_speed=kps4.set.v_reel_out, dt)
end
lift, drag = lift_drag(kps4)
sys_state = KiteModels.SysState(kps4)
log!(logger, sys_state)
elev = rad2deg(logger.elevation_vec[end])
println("Lift: $lift, Drag: $drag, elev: $elev, Iterations: $(kps4.iter)")

q = QuatRotation(calc_orient_quat(kps4; viewer=false))
roll, pitch, yaw = rad2deg.(quat2euler(q))
println("--> orient_quat:       roll: ", roll, " pitch: ", pitch, "  yaw: ", yaw)
roll, pitch, yaw = rad2deg.(orient_euler(kps4))
println("--> orient_euler:      roll: ", roll, " pitch: ", pitch, " yaw:  ", yaw)
q = QuatRotation(calc_orient_quat(kps4; viewer=true))
roll, pitch, yaw = rad2deg.(quat2euler(q))
println("--> orient_quat (viewer): roll: ", roll, " pitch: ", pitch, "   yaw: ", yaw)

println("x:", kps4.x) # from trailing edge to leading edge in ENU reference frame
println("y:", kps4.y) # to the right looking in flight direction
println("z:", kps4.z) # down
azimuth = calc_azimuth(kps4)
println("azimuth: ", round(rad2deg(azimuth), digits = 2), "°")
println("azimuth_north: ", round(rad2deg(KiteUtils.azimuth_north(pos_kite(kps4))), digits = 2), "°")

# print point C and point D
pos_C, pos_D = kps4.pos[kps4.set.segments+4], kps4.pos[kps4.set.segments+5]

# print alpha2, alpha3, alpha4
println("alpha2, alpha3, alpha4: ", kps4.alpha_2, " ", kps4.alpha_3, " ", kps4.alpha_4)
println("heading: ", round(rad2deg(calc_heading(kps4)), digits = 2), "°")

ss = SysState(kps4)
println("AoA: ", rad2deg(ss.AoA), "°")
println("CL2: ", ss.CL2, " CD2: ", ss.CD2)
println("v_wind_gnd: ", ss.v_wind_gnd)
println("v_wind_200m: ", ss.v_wind_200m)
println("v_wind_kite: ", ss.v_wind_kite)

# output on main branch
# Lift: 1047.1795339611076, Drag: 281.39765463928745, elev: 72.77014, Iterations: 1552
# roll: -0.0 pitch: 5.755691707832273 yaw: 90.00000250447816
# x:[-0.9949585608372032, 0.0, 0.10028689952711267]
# y:[0.0, 1.0, 0.0]
# z:[-0.10028689952711267, -2.8981421002141166e-19, -0.9949585608372032]
# alpha2, alpha3, alpha4: 7.546125780476343 10.000000853773646 10.000000853773646

