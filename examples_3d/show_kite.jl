# generate test cases for the calculation of roll, pitch and yaw
using LinearAlgebra
import ReferenceFrameRotations as RFR

using Printf

using Pkg, Timers
tic()
if ! ("KiteViewers" ∈ keys(Pkg.project().dependencies))
    Pkg.activate("examples_3d")
    pkg"add KiteModels#main"
end
using KiteUtils, Rotations, StaticArrays
using KiteViewers
toc()

yaw = deg2rad(0)
pitch = deg2rad(0)
roll = deg2rad(0)

viewer::Viewer3D = Viewer3D(true);
segments=6
state=demo_state_4p(segments+1, 6)

"""
    is_right_handed_orthonormal(x, y, z)

Returns `true` if the vectors `x`, `y` and `z` form a right-handed orthonormal basis.
"""
function is_right_handed_orthonormal(x, y, z)
    R = [x y z]
    R*R' ≈ I && det(R) ≈ 1
end

"""
    rot3d(ax, ay, az, bx, by, bz)

Calculate the rotation matrix that needs to be applied on the reference frame (ax, ay, az) to match 
the reference frame (bx, by, bz).
All parameters must be 3-element vectors. Both refrence frames must be orthogonal,
all vectors must already be normalized.

Source: [TRIAD_Algorithm](http://en.wikipedia.org/wiki/User:Snietfeld/TRIAD_Algorithm)
"""
function rot3d(ax, ay, az, bx, by, bz)
    @assert is_right_handed_orthonormal(ax, ay, az)
    @assert is_right_handed_orthonormal(bx, by, bz)    
    R_ai = [ax az ay]
    R_bi = [bx bz by]
    return R_bi * R_ai'
end

function calc_orient_rot(x, y, z; old=false)
    if old
        pos_kite_ = ones(3)
        pos_before = pos_kite_ .+ z
    
        rotation = rot(pos_kite_, pos_before, -x)
    else
        # reference frame for the orientation: NED (north, east, down)
        ax = [0, 1, 0] # in ENU reference frame this is pointing to the south
        ay = [1, 0, 0] # in ENU reference frame this is pointing to the west
        az = [0, 0, -1] # in ENU reference frame this is pointing down
        rotation = rot3d(ax, ay, az, x, y, z)
    end
    return rotation
end

quat2frame(q::AbstractMatrix) = quat2frame(QuatRotation(q))
function quat2frame(q::QuatRotation)
    x = [0,  1.0, 0]
    y = [1.0,  0, 0]
    z = [0,    0, -1.0]
    return q*x, q*y, q*z
end

function new2old(q::QuatRotation)
    x, y, z = quat2frame(q)
    rot = calc_orient_rot(x, y, z; old=true)
    q = QuatRotation(rot)
    return Rotations.params(q)
end

# z-y′-x″ (intrinsic rotations) or x-y-z (extrinsic rotations): 
# the intrinsic rotations are known as: yaw, pitch and roll

# x: from trailing edge to leading edge
# y: to the right looking in flight direction
# z: down

# If kite (x axis) is pointing to the north, and is at zenith, then in ENUs reference frame:
# - x = 0, 1, 0
# - y = 1, 0, 0
# - z = 0, 0,-1
# This would be equal to the NED reference frame.
x = [0, 1, 0]
y = [1, 0, 0]
z = [0, 0,-1]

# R = Yaw * Pitch * Roll

D1 = RFR.angle_to_dcm(yaw, pitch, roll, :ZYX)
x4 = D1 * x
y4 = D1 * y
z4 = D1 * z
println("Yaw: ", rad2deg(yaw), ", Pitch: ", rad2deg(pitch), ", Roll: ", rad2deg(roll), 
        "\nx = ", x4, "\ny = ", y4, "\nz = ", z4)

euler = RFR.dcm_to_angle(D1, :ZYX)
yaw = euler.a1
pitch = euler.a2
roll = euler.a3
println("Yaw: ", rad2deg(yaw), ", Pitch: ", rad2deg(pitch), ", Roll: ", rad2deg(roll))

q = QuatRotation(calc_orient_rot(x4, y4, z4))
roll, pitch, yaw = rad2deg.(quat2euler(q))

state.orient .= new2old(q)
update_system(viewer, state, kite_scale=0.25)