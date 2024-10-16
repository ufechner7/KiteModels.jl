# generate test cases for the calculation of roll, pitch and yaw
using LinearAlgebra
import ReferenceFrameRotations as RFR

using Printf

using Pkg, Timers
tic()
if ! ("KiteViewers" ∈ keys(Pkg.project().dependencies))
    Pkg.activate("examples_3d")
    pkg"add KiteUtils#main"
    pkg"add KiteModels#main"
end
using KiteUtils, Rotations, StaticArrays
using KiteViewers
toc()

# yaw = deg2rad(-63.529095)
# pitch = deg2rad(9.046745)
# roll = deg2rad(3.800396)
# yaw = deg2rad(0)   # noise pointing to the north
# yaw = deg2rad(180) # noise pointing to the south
yaw = deg2rad(-90)  # noise pointing to the west
pitch = deg2rad(0)
roll = deg2rad(0)


# x: from trailing edge to leading edge
# y: to the right looking in flight direction
# z: down

# reference frame in NED
x = [1,  0, 0] # nose pointing to the north
y = [0,  1, 0] # wing pointing to the east
z = [0,  0, 1] # z pointing down

"""
    is_right_handed_orthonormal(x, y, z)

Returns `true` if the vectors `x`, `y` and `z` form a right-handed orthonormal basis.
"""
function is_right_handed_orthonormal(x, y, z)
    R = [x y z]
    R*R' ≈ I && det(R) ≈ 1
end

function euler2rot(roll, pitch, yaw)
    φ      = roll
    R_x = [1    0       0;
              0  cos(φ) -sin(φ);
              0  sin(φ)  cos(φ)]
    θ      = pitch          
    R_y = [ cos(θ)  0  sin(θ);
                 0     1     0;
              -sin(θ)  0  cos(θ)]
    ψ      = yaw
    R_z = [cos(ψ) -sin(ψ) 0;
              sin(ψ)  cos(ψ) 0;
                 0       0   1]
    R   = R_z * R_y * R_x
    return R
end

"""
    quat2viewer(q::QuatRotation)
    quat2viewer(rot::AbstractMatrix)
    quat2viewer(orient::AbstractVector)

Convert the quaternion q to the viewer reference frame. It can also be passed
as a rotation matrix or as 4-element vector [w,i,j,k], where w is the real part
and i, j, k are the imaginary parts of the quaternion.
"""
quat2viewer_(rot::AbstractMatrix) = quat2viewer(QuatRotation(rot))
quat2viewer_(orient::AbstractVector) = quat2viewer(QuatRotation(orient))
function quat2viewer_(q::QuatRotation)
    # 1. get reference frame
    rot = inv(RotMatrix{3}(q))
    x = enu2ned(rot[1,:])
    y = enu2ned(rot[2,:])
    z = enu2ned(rot[3,:])
    # 2. convert it using the old method
    ax = [0, 1, 0] # in ENU reference frame this is pointing to the south
    ay = [1, 0, 0] # in ENU reference frame this is pointing to the west
    az = [0, 0, -1] # in ENU reference frame this is pointing down
    rotation = rot3d(ax, ay, az, x, y, z)
    q_old = QuatRotation(rotation)
    x = [0,  1.0, 0]
    y = [1.0,  0, 0]
    z = [0,    0, -1.0]
    x, y, z = q_old*x, q_old*y, q_old*z
    rot = calc_orient_rot(x, y, z; viewer=true, ENU=false)
    q = QuatRotation(rot)
    return Rotations.params(q)
end

D1 = euler2rot(roll, pitch, yaw)
x4 = D1 * x
y4 = D1 * y
z4 = D1 * z
println("Yaw: ", rad2deg(yaw), ", Pitch: ", rad2deg(pitch), ", Roll: ", rad2deg(roll), 
        "\nNED:\nx4 = ", (x4), "\ny4 = ", (y4), "\nz4 = ", (z4))

roll, pitch, yaw = quat2euler(QuatRotation(D1))
println("Yaw: ", rad2deg(yaw), ", Pitch: ", rad2deg(pitch), ", Roll: ", rad2deg(roll))

rot = calc_orient_rot(x4, y4, z4; ENU=false)
q = QuatRotation(rot)

viewer::Viewer3D = Viewer3D(true);
segments=6
state=demo_state_4p(segments+1, 12; yaw=pi-yaw)
state.orient = quat2viewer_(q)
update_system(viewer, state, kite_scale=0.25)
nothing