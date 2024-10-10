# test calculation of the orientation, kite pointing to the west and is at zenith
import KiteUtils
using LinearAlgebra, Rotations

# If kite (x-axis) is pointing to the west, against the wind, and is at zenith, then:
# - x = -1, 0, 0
# - y =  0, 1, 0
# - z =  0, 0,-1

x = [-1, 0, 0]
y = [ 0, 1, 0]
z = [ 0, 0,-1]

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
roll, pitch, yaw = rad2deg.(KiteUtils.quat2euler(q))
println("--> orient_quat:       roll: ", roll, " pitch:  ", pitch, "  yaw: ", yaw)

q = Rotations.params(q)
roll, pitch, yaw = rad2deg.(KiteUtils.quat2euler(q))
println("--> orient_quat:       roll: ", roll, " pitch:  ", pitch, "  yaw: ", yaw)
rot

