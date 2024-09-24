# plot the lift and drag coefficients as function of angle of attack

using Printf
using KiteModels, KitePodModels, KiteUtils, LinearAlgebra, Rotations

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

function quat2euler(q)
    # Convert quaternion to RotXYZ
    rot = RotXYZ(q)
    
    # Extract roll, pitch, and yaw from RotXYZ
    roll = rot.theta1
    pitch = rot.theta2
    yaw = rot.theta3

    return roll, pitch, yaw
end

elev = set.elevation
i = 1
set.v_wind = V_WIND # 25
logger = Logger(set.segments + 5, STEPS)

kcu::KCU = KCU(set)
kps4::KPS4 = KPS4(kcu)
integrator = KiteModels.init_sim!(kps4; delta=0.03, stiffness_factor=0.1, prn=STATISTIC)
lift, drag = lift_drag(kps4)
sys_state = KiteModels.SysState(kps4)
log!(logger, sys_state)
elev = rad2deg(logger.elevation_vec[end])
println("Lift: $lift, Drag: $drag, elev: $elev, Iterations: $(kps4.iter)")

q = QuatRotation(sys_state.orient)
# println(q)
roll, pitch, yaw = rad2deg.(quat2euler(q))
println("roll: ", roll, " pitch: ", pitch, " yaw: ", yaw)
println("x:", kps4.x)
println("y:", kps4.y)
println("z:", kps4.z)

# print point C and point D
pos_C, pos_D = kps4.pos[kps4.set.segments+4], kps4.pos[kps4.set.segments+5]

# output on main branch
# Lift: 1047.1795339611076, Drag: 281.39765463928745, elev: 72.77014, Iterations: 1552
# roll: -0.0 pitch: 5.755691707832273 yaw: 90.00000250447816
# x:[-0.9949585608372032, 0.0, 0.10028689952711267]
# y:[0.0, 1.0, 0.0]
# z:[-0.10028689952711267, -2.8981421002141166e-19, -0.9949585608372032]
