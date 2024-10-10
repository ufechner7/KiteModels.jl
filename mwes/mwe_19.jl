# test calculation of the orientation
using LinearAlgebra, Rotations

# z-y′-x″ (intrinsic rotations) or x-y-z (extrinsic rotations): 
# the intrinsic rotations are known as: yaw, pitch and roll

# x: from trailing edge to leading edge
# y: to the right looking in flight direction
# z: down

# If x, y and z are given in ENU
# x = [0, 1, 0] y = [0, 0, 1] z = [1, 0, 0] should give 90 degrees pitch
x = [ 0, 0, 1]
y = [ 1, 0, 0]
z = [ 0, 1, 0]

@assert is_right_handed_orthonormal(x, y, z)

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

quat2euler(q::AbstractVector) = quat2euler(QuatRotation(q))
function quat2euler(q::QuatRotation)
    # Convert quaternion to RotXYZ
    rot = RotXYZ(q)
    
    # Extract roll, pitch, and yaw from RotXYZ
    roll = rot.theta1
    pitch = rot.theta2
    yaw = rot.theta3

    return roll, pitch, yaw
end

function calc_orient_quat(x, y, z)
    # reference frame for the orientation: NED
    ax = [0, 1,  0] # in ENU reference frame this is pointing to the north
    ay = [1, 0,  0] # in ENU reference frame this is pointing to the east
    az = [0, 0, -1] # in ENU reference frame this is pointing down
    @assert is_right_handed_orthonormal(ax, ay, az)
    rot = rot3d(ax, ay, az, x, y, z)
    return rot
    # q = QuatRotation(rotation)
    # return Rotations.params(q)
end
rot = calc_orient_quat(x, y, z)
q = QuatRotation(rot)
println("q: ", q)
roll, pitch, yaw = rad2deg.(quat2euler(q))
println("--> orient_quat:       roll: ", roll, " pitch:  ", pitch, "  yaw: ", yaw)
@test pitch ≈ 90

q = Rotations.params(q)
roll, pitch, yaw = rad2deg.(quat2euler(q))
println("--> orient_quat:       roll: ", roll, " pitch:  ", pitch, "  yaw: ", yaw)
rot

