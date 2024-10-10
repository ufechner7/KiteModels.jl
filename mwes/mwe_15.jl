# test calculation of the orientation
import KiteUtils
using LinearAlgebra, Rotations

# If kite (x axis) is pointing to the north, and is at zenith, then:
# - x = 0, 1, 0
# - y = 1, 0, 0
# - z = 0, 0,-1
# This would be the NED reference frame.

x = [0, 1, 0]
y = [1, 0, 0]
z = [0, 0,-1]

function calc_orient_quat(x, y, z)
    # reference frame for the orientation: SWD (south, west, down)
    ax = [0, -1, 0] # in ENU reference frame this is pointing to the south
    ay = [-1, 0, 0] # in ENU reference frame this is pointing to the west
    az = [0, 0, -1] # in ENU reference frame this is pointing down
    rotation = KiteUtils.rot3d(ax, ay, az, x, y, z)

    q = QuatRotation(rotation)
    return Rotations.params(q)
end

q = calc_orient_quat(x, y, z)
roll, pitch, yaw = rad2deg.(KiteUtils.quat2euler(q))
println("--> orient_quat:       roll: ", roll, " pitch:  ", pitch, "  yaw: ", yaw)

