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
    R_ai = hcat(ax, az, ay)
    R_bi = hcat(bx, bz, by)
    return R_bi * R_ai'
end

ax = [1, 0, 0] 
ay = [0, 1, 0] 
az = [0, 0, 1] 

x = [1, 0, 0] 
y = [0, 1, 0] 
z = [0, 0, 1]

rot1 = rot3d(ax, ay, az, x, y, z)
q1 = QuatRotation(rot1)
@test all(Rotations.params(q) .== SVector{4, Float64}([1.0 0 0 0]))

x = [-1, 0, 0] 
y = [0, 1, 0] 
z = [0, 0, 1]
rot2 = rot3d(ax, ay, az, x, y, z)
q2 = QuatRotation(rot2)
