# test calculation of the orientation, kite pointing to the west and is at zenith
import KiteUtils
using LinearAlgebra, Rotations

# z-y′-x″ (intrinsic rotations) or x-y-z (extrinsic rotations): 
# the intrinsic rotations are known as: yaw, pitch and roll

# If kite (x-axis) is pointing to the west, against the wind, and is at zenith, then:
# - x = -1, 0, 0
# - y =  0, 1, 0
# - z =  0, 0,-1

x = [-1, 0, 0]
y = [ 0, 1, 0]
z = [ 0, 0,-1]
quat2euler(q::AbstractVector) = quat2euler(QuatRotation(q))
function quat2euler(q::QuatRotation)
    # Convert quaternion to RotXYZ
    rot = RotXYZ(q)
    
    # Extract roll, pitch, and yaw from RotXYZ
    roll = rot.theta1
    pitch = rot.theta2
    yaw = rot.theta3

    return roll, pitch, yaw
end

function calc_orient_quat(x, y, z)
    # reference frame for the orientation: NED
    ax = [0, 1,  0] # in ENU reference frame this is pointing to the north
    ay = [1, 0,  0] # in ENU reference frame this is pointing to the east
    az = [0, 0, -1] # in ENU reference frame this is pointing down
    rot = KiteUtils.rot3d(ax, ay, az, x, y, z)
    return rot
    # q = QuatRotation(rotation)
    # return Rotations.params(q)
end
rot = calc_orient_quat(x, y, z)
q = QuatRotation(rot)
roll, pitch, yaw = rad2deg.(quat2euler(q))
println("--> orient_quat:       roll: ", roll, " pitch:  ", pitch, "  yaw: ", yaw)

q = Rotations.params(q)
roll, pitch, yaw = rad2deg.(quat2euler(q))
println("--> orient_quat:       roll: ", roll, " pitch:  ", pitch, "  yaw: ", yaw)
rot

