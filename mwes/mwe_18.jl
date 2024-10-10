# test calculation of the orientation
using LinearAlgebra, Rotations, StaticArrays, Test
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
    R_ai = hcat(ax, az, ay)
    R_bi = hcat(bx, bz, by)
    return R_bi * R_ai'
end

function is_right_handed_orthonormal(x, y, z)
    if !(norm(x) ≈ 1) || !(norm(y) ≈ 1) || !(norm(z) ≈ 1)
        return false
    end
    if !((x ⋅ y) ≈ 0) || !((y ⋅ z) ≈ 0) || !((z ⋅ x) ≈ 0)
        return false
    end
    return det([x y z]) ≈ 1
end

ax = [1, 0, 0] 
ay = [0, 1, 0] 
az = [0, 0, 1]
@assert is_right_handed_orthonormal(ax, ay, az)

x = [0, 1, 0]
y = [1, 0, 0]
z = [0, 0,-1]
@assert is_right_handed_orthonormal(x, y, z)

rot1 = rot3d(ax, ay, az, x, y, z)
q1 = QuatRotation(rot1)

