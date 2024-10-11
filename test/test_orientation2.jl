# unit tests for calculation of the orientation
using LinearAlgebra, Rotations, Test, StaticArrays

# Kite reference frame
# x: from trailing edge to leading edge
# y: to the right looking in flight direction
# z: down

# all coordinates are in ENU (East, North, Up) reference frame
# the orientation is calculated with respect to the NED (North, East, Down) reference frame

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

function calc_orient_rot(x, y, z)
    # reference frame for the orientation: NED
    ax = @SVector[0, 1, 0] # in ENU reference frame this is pointing to the north
    ay = @SVector[1, 0, 0] # in ENU reference frame this is pointing to the east
    az = @SVector[0, 0,-1] # in ENU reference frame this is pointing down
    rot = rot3d(ax, ay, az, x, y, z)
    return rot
end

quat2euler(q::AbstractVector) = quat2euler(QuatRotation(q))
function quat2euler(q::QuatRotation)
    # Convert quaternion to RotXYZ
    rot = RotXYZ(q)
    
    # Extract roll, pitch, and yaw from RotXYZ
    roll = rot.theta2
    pitch = rot.theta1
    yaw = -rot.theta3

    return roll, pitch, yaw
end

@testset "calc_orientation, kite pointing to the north and is at zenith" begin
    # If kite (x axis) is pointing to the north, and is at zenith, then in ENUs reference frame:
    # - x = 0, 1, 0
    # - y = 1, 0, 0
    # - z = 0, 0,-1
    # This would be equal to the NED reference frame.
    x = [0, 1, 0]
    y = [1, 0, 0]
    z = [0, 0,-1]
    @assert is_right_handed_orthonormal(x, y, z)
    rot = calc_orient_rot(x, y, z)
    q = QuatRotation(rot)
    roll, pitch, yaw = rad2deg.(quat2euler(q))
    @test roll ≈ 0
    @test pitch ≈ 0
    @test yaw ≈ 0
end
@testset "calc_orientation, kite pointing to the west and is at zenith " begin
    x = [-1, 0, 0]
    y = [ 0, 1, 0]
    z = [ 0, 0,-1]
    @assert is_right_handed_orthonormal(x, y, z)
    rot = calc_orient_rot(x, y, z)
    q = QuatRotation(rot)
    roll, pitch, yaw = rad2deg.(quat2euler(q))
    @test roll ≈ 0
    @test pitch ≈ 0
    @test yaw ≈ -90
end
@testset "calc_orientation, kite pointing to the north, right tip up   " begin
    # x = [0, 1, 0] y = [0, 0, 1] z = [1, 0, 0] should give -90 degrees roll
    x = [ 0, 1, 0]
    y = [ 0, 0, 1]
    z = [ 1, 0, 0]
    @assert is_right_handed_orthonormal(x, y, z)
    rot = calc_orient_rot(x, y, z)
    q = QuatRotation(rot)
    roll, pitch, yaw = rad2deg.(quat2euler(q))
    @test roll ≈ -90
    @test pitch ≈ 0
    @test yaw ≈ 0
end
@testset "calc_orientation, kite pointing upwards, right tip eastwards " begin
    # If x, y and z are given in ENU
    # x = [0, 0, 1] y = [1, 0, 0] z = [0, 1, 0] should give 90 degrees pitch
    x = [ 0, 0, 1]  # nose pointing up
    y = [ 1, 0, 0]  # right tip pointing east
    z = [ 0, 1, 0]  # z axis pointing to the north
    @assert is_right_handed_orthonormal(x, y, z)
    rot = calc_orient_rot(x, y, z)
    q = QuatRotation(rot)
    roll, pitch, yaw = rad2deg.(quat2euler(q))
    @test roll ≈ 0
    @test pitch ≈ 90
    @test yaw ≈ 0
end
# @testset "calc_orientation, all angles positive                        " begin
#     # x, y and z are given in ENU
#     x = [ 0.61237244,  0.61237244,  0.5       ]
#     y = [ 0.65973961, -0.04736717, -0.75      ]
#     z = [-0.43559574,  0.78914913, -0.4330127 ]
#     @assert is_right_handed_orthonormal(x, y, z)
#     rot = calc_orient_rot(x, y, z)
#     q = QuatRotation(rot)
#     roll, pitch, yaw = rad2deg.(quat2euler(q))
#     @test_broken roll ≈ 60
#     @test_broken pitch ≈ 30
#     @test_broken yaw ≈ 45
#     println("roll: $roll, pitch: $pitch, yaw: $yaw")
# end
@testset "calc_orientation, yaw = 20°                                  " begin
    # x, y and z are given in ENU
    x =  [0.34202014332566866, 0.9396926207859083, 0.0]
    y = [0.9396926207859083, -0.34202014332566866, 0.0]
    z = [ 0,  0,  -1 ]
    @assert is_right_handed_orthonormal(x, y, z)
    rot = calc_orient_rot(x, y, z)
    q = QuatRotation(rot)
    roll, pitch, yaw = rad2deg.(quat2euler(q))
    @test roll  ≈ 0
    @test pitch ≈ 0
    @test yaw   ≈ 20
end
@testset "calc_orientation, pitch = 30°                                " begin
    global rot
    # x, y and z are given in ENU
    # x, y and z are given in ENU
    x = [0.0, 0.8660254037844387, 0.49999999999999994]
    y = [ 1.,         0.,                 0.       ]
    z = [0.0, 0.49999999999999994, -0.8660254037844387]
    @assert is_right_handed_orthonormal(x, y, z)
    rot = calc_orient_rot(x, y, z)
    q = QuatRotation(rot)
    roll, pitch, yaw = rad2deg.(quat2euler(q))
    println("roll: $roll, pitch: $pitch, yaw: $yaw")
    @test roll  ≈ 0
    @test pitch ≈ 30.0
    @test yaw   ≈ 0
end
nothing
