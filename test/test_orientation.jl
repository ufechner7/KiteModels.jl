# unit tests for calculation of the orientation
using LinearAlgebra, Rotations, Test, StaticArrays
import ReferenceFrameRotations as RFR
using Pkg
pkg"add KiteUtils#main"
using KiteUtils

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

function enu2ned(vec::AbstractVector)  
    R = [0 1 0; 1 0 0; 0 0 -1]
    R*vec
end

quat2euler2(q::AbstractVector) = quat2euler(QuatRotation(q))
function quat2euler2(q::QuatRotation)  
    D = RFR.DCM(q)
    pitch = asin(−D[3,1])
    roll  = atan(D[3,2], D[3,3])
    yaw   = atan(D[2,1], D[1,1])
    return roll, pitch, yaw
end

function calc_orient_rot2(x, y, z; viewer=false)
    if viewer
        pos_kite_ = @SVector ones(3)
        pos_before = pos_kite_ .+ z
        rotation = rot(pos_kite_, pos_before, -x)
    else
        # reference frame for the orientation: NED (north, east, down)
        ax = @SVector [1, 0, 0]
        ay = @SVector [0, 1, 0]
        az = @SVector [0, 0, 1]
        rotation = rot3d(ax, ay, az, x, y, z)
    end
    return rotation
end

@testset "calc_orientation, kite pointing to the north and is at zenith" begin
    # If kite (x axis) is pointing to the north, and is at zenith, then in ENUs reference frame:
    # - x = 0, 1, 0
    # - y = 1, 0, 0
    # - z = 0, 0,-1
    # This would be equal to the NED reference frame.
    x = enu2ned([0, 1, 0])
    y = enu2ned([1, 0, 0])
    z = enu2ned([0, 0,-1])
    @assert is_right_handed_orthonormal(x, y, z)
    rot = calc_orient_rot2(x, y, z)
    q = QuatRotation(rot)
    roll, pitch, yaw = rad2deg.(quat2euler2(q))
    @test roll ≈ 0
    @test pitch ≈ 0
    @test yaw ≈ 0
end
@testset "calc_orientation, kite pointing to the west and is at zenith " begin
    x = [0.0, -1.0, 0.0] # NED, kite pointing to the west
    y = [1.0, 0.0, 0.0]  # NED, right tip pointing to the north
    z = [0.0, 0.0, 1.0]  # NED, z axis pointing down
    @assert is_right_handed_orthonormal(x, y, z)
    rot = calc_orient_rot2(x, y, z)
    q = QuatRotation(rot)
    roll, pitch, yaw = rad2deg.(quat2euler2(q))
    @test roll ≈ 0
    @test pitch ≈ 0
    @test yaw ≈ -90
end
@testset "calc_orientation, kite pointing to the north, right tip up   " begin
    # x = [0, 1, 0] y = [0, 0, 1] z = [1, 0, 0] should give -90 degrees roll
    x = enu2ned([ 0, 1, 0])
    y = enu2ned([ 0, 0, 1])
    z = enu2ned([ 1, 0, 0])
    @assert is_right_handed_orthonormal(x, y, z)
    rot = calc_orient_rot2(x, y, z)
    q = QuatRotation(rot)
    roll, pitch, yaw = rad2deg.(quat2euler2(q))
    @test roll ≈ -90
    @test pitch ≈ 0
    @test yaw ≈ 0
end
@testset "calc_orientation, kite pointing upwards, right tip eastwards " begin
    # If x, y and z are given in ENU
    # x = [0, 0, 1] y = [1, 0, 0] z = [0, 1, 0] should give 90 degrees pitch
    x = enu2ned([ 0, 0, 1])  # nose pointing up
    y = enu2ned([ 1, 0, 0])  # right tip pointing east
    z = enu2ned([ 0, 1, 0])  # z axis pointing to the north
    @assert is_right_handed_orthonormal(x, y, z)
    rot = calc_orient_rot2(x, y, z)
    q = QuatRotation(rot)
    roll, pitch, yaw = rad2deg.(quat2euler2(q))
    @test roll ≈ 0
    @test pitch ≈ 90
    @test yaw ≈ 0
end
@testset "calc_orientation, all angles positive                        " begin
    # x, y and z are given in ENU
    x = [0.8137976813493738, 0.29619813272602386, -0.49999999999999994]
    y = [0.0400087565481419, 0.8297694655894313, 0.5566703992264194]
    z = [0.5797694655894312, -0.47302145844036114, 0.6634139481689384]
    @assert is_right_handed_orthonormal(x, y, z)
    rot = calc_orient_rot2(x, y, z)
    q = QuatRotation(rot)
    roll, pitch, yaw = rad2deg.(quat2euler2(q))
    @test roll ≈ 40
    @test pitch ≈ 30
    @test yaw ≈ 20
end
# @testset "calc_orientation, yaw = 20°                                  " begin
#     # x, y and z are given in ENU
#     x = enu2ned([0.34202014332566866, 0.9396926207859083, 0.0])
#     y = enu2ned([0.9396926207859083, -0.34202014332566866, 0.0])
#     z = enu2ned([ 0,  0,  -1 ])
#     @assert is_right_handed_orthonormal(x, y, z)
#     rot = calc_orient_rot2(x, y, z)
#     q = QuatRotation(rot)
#     roll, pitch, yaw = rad2deg.(quat2euler2(q))
#     @test roll  ≈ 0
#     @test pitch ≈ 0
#     @test yaw   ≈ 20
# end
# @testset "calc_orientation, pitch = 30°                                " begin
#     # x, y and z are given in ENU
#     x = enu2ned([0.0, 1.0, 0.0])
#     y = enu2ned([0.8660254037844385, 0.0, 0.5000000000000002])
#     z = enu2ned([0.5000000000000002, 0.0, -0.8660254037844385])
#     @assert is_right_handed_orthonormal(x, y, z)
#     rot = calc_orient_rot2(x, y, z)
#     q = QuatRotation(rot)
#     roll, pitch, yaw = rad2deg.(quat2euler2(q))
#     @test roll  ≈ 0
#     @test pitch ≈ 30.0
#     @test yaw   ≈ 0
# end
# @testset "calc_orientation, yaw=20°, pitch = 30°                       " begin
#     # x, y and z are given in ENU
#     x = [0.29619813272602386, 0.9396926207859083, 0.1710100716628345]
#     y = [0.8137976813493735, -0.3420201433256688, 0.4698463103929544]
#     z = [0.5000000000000002, 1.2018516789897272e-17, -0.8660254037844385]
    
#     @assert is_right_handed_orthonormal(x, y, z)
#     rot = calc_orient_rot(x, y, z)
#     q = QuatRotation(rot)
#     roll, pitch, yaw = rad2deg.(quat2euler(q))
#     @test roll + 1  ≈ 0 + 1
#     @test pitch     ≈ 30.0
#     @test yaw       ≈ 20.0
# end
# @testset "calc_orientation, roll = 40°                                 " begin
#     # x, y and z are given in ENU
#     x = [0.0, 0.7660444431189782, -0.6427876096865391]
#     y = [1.0, 0.0, 0.0]
#     z = [0.0, -0.6427876096865391, -0.7660444431189782]
    
#     @assert is_right_handed_orthonormal(x, y, z)
#     rot = calc_orient_rot(x, y, z)
#     q = QuatRotation(rot)
#     roll, pitch, yaw = rad2deg.(quat2euler(q))
#     @test roll       ≈ 40
#     @test pitch+1    ≈ 0.0+1
#     @test yaw+1      ≈ 0.0+1
# end
nothing
