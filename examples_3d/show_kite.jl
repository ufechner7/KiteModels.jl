# Copyright (c) 2024 Uwe Fechner
# SPDX-License-Identifier: MIT

# show a static kite for verification of roll, pitch and yaw
# To run this script, use the following commands:
# julia --project
# include("examples_3d/show_kite.jl")
using LinearAlgebra

using Printf

using Pkg, Timers
tic()
if ! ("KiteViewers" ∈ keys(Pkg.project().dependencies))
    Pkg.activate("examples_3d")
    # pkg"add KiteUtils#main"
    # pkg"add KiteModels#main"
end
using KiteUtils, Rotations, StaticArrays
using KiteViewers
toc()

# yaw = deg2rad(0)   # noise pointing to the north
yaw = deg2rad(180) # noise pointing to the south
# yaw = deg2rad(-90)  # noise pointing to the west
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

D1 = euler2rot(roll, pitch, yaw)
x4 = D1 * x
y4 = D1 * y
z4 = D1 * z
println("Yaw: ", rad2deg(yaw), ", Pitch: ", rad2deg(pitch), ", Roll: ", rad2deg(roll), 
        "\nNED:\nx4 = ", (x4), "\ny4 = ", (y4), "\nz4 = ", (z4))

roll, pitch, yaw = quat2euler(QuatRotation(D1))
println("Yaw: ", rad2deg(yaw), ", Pitch: ", rad2deg(pitch), ", Roll: ", rad2deg(roll))
println("azimuth_north: ", rad2deg(pi-yaw))

rot = calc_orient_rot(x4, y4, z4; ENU=false)
q = QuatRotation(rot)

viewer::Viewer3D = Viewer3D(true);
segments=6
state=demo_state_4p(segments+1, 12; azimuth_north=pi-yaw)
state.orient = quat2viewer(q)
update_system(viewer, state, kite_scale=0.25)
nothing