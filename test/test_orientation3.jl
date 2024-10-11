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
    rot = rot3d(x, y, z, ax, ay, az)
    return rot
end

@testset "calc_orientation, pitch = 30°                                " begin
    global rot
    # x, y and z are given in ENU
    x = [0.0, 0.8660254037844387, 0.49999999999999994]
    y = [ 1.,         0.,                 0.       ]
    z = [0.0, 0.49999999999999994, -0.8660254037844387]
    @assert is_right_handed_orthonormal(x, y, z)
    rot = calc_orient_rot(x, y, z)
end
rot
