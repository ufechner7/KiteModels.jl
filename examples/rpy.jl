using LinearAlgebra, Rotations, Test, StaticArrays
import ReferenceFrameRotations as RFR
using Pkg
pkg"add KiteUtils#main"
using KiteUtils


"""
    is_right_handed_orthonormal(x, y, z)

Returns `true` if the vectors `x`, `y` and `z` form a right-handed orthonormal basis.
"""
function is_right_handed_orthonormal(x, y, z)
    R = [x y z]
    R*R' ≈ I && det(R) ≈ 1
end

function xyz_to_euler(x_, y_, z_)
    x = -x_
    y = y_
    z = -z_
    rot = rot3d([1, 0, 0], [0, 1, 0], [0, 0, 1], x, y, z)
    q = QuatRotation(rot).q

    # ref: https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
    sinr_cosp = 2 * (q.s * q.v1 + q.v2 * q.v3)
    cosr_cosp = 1 - 2 * (q.v1 * q.v1 + q.v2 * q.v2)
    roll = atan(sinr_cosp, cosr_cosp)

    # pitch (y-axis rotation)
    sinp = sqrt(1 + 2 * (q.s * q.v2 - q.v1 * q.v3))
    cosp = sqrt(1 - 2 * (q.s * q.v2 - q.v1 * q.v3))
    pitch = 2 * atan(sinp, cosp) - π / 2

    # yaw (z-axis rotation)
    siny_cosp = 2 * (q.s * q.v3 + q.v1 * q.v2)
    cosy_cosp = 1 - 2 * (q.v2 * q.v2 + q.v3 * q.v3)
    yaw = atan(siny_cosp, cosy_cosp)

    return roll, pitch, yaw
end

# nothing
x = [ -1, 0, 0]
y = [ 0, 1, 0]
z = [ 0, 0, -1]
@show rad2deg.(xyz_to_euler(x, y, z))

println("roll 90 deg")
x = [ -1, 0, 0]
y = [ 0, 0, 1]
z = [ 0, 1, 0]
@show rad2deg.(xyz_to_euler(x, y, z))

println("pitch 90 deg")
x = [ 0, 0, 1]
y = [ 0, 1, 0]
z = [ -1, 0, 0]
@show rad2deg.(xyz_to_euler(x, y, z))

println("yaw 90 deg")
x = [ 0, -1, 0]
y = [ -1, 0, 0]
z = [ 0, 0, -1]
@show rad2deg.(xyz_to_euler(x, y, z))
nothing