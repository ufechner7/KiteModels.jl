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

yaw = deg2rad(-90)
pitch = deg2rad(0)
roll = deg2rad(0)

x = enu2ned([0, 1, 0])
y = enu2ned([1, 0, 0])
z = enu2ned([0, 0,-1])

"""
    is_right_handed_orthonormal(x, y, z)

Returns `true` if the vectors `x`, `y` and `z` form a right-handed orthonormal basis.
"""
function is_right_handed_orthonormal(x, y, z)
    R = [x y z]
    R*R' ≈ I && det(R) ≈ 1
end

# x: from trailing edge to leading edge
# y: to the right looking in flight direction
# z: down

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

D1 = euler2rot(roll, pitch, yaw)
x4 = D1 * x
y4 = D1 * y
z4 = D1 * z
println("Yaw: ", rad2deg(yaw), ", Pitch: ", rad2deg(pitch), ", Roll: ", rad2deg(roll), 
        "\nx = ", ned2enu(x4), "\ny = ", ned2enu(y4), "\nz = ", ned2enu(z4))

roll, pitch, yaw = quat2euler(QuatRotation(D1))
println("Yaw: ", rad2deg(yaw), ", Pitch: ", rad2deg(pitch), ", Roll: ", rad2deg(roll))

rot = calc_orient_rot(x4, y4, z4; ENU=true)
q = QuatRotation(rot)

viewer::Viewer3D = Viewer3D(true);
segments=6
state=demo_state_4p(segments+1, 12; yaw)
state.orient .= quat2viewer(q)
update_system(viewer, state, kite_scale=0.25)
nothing