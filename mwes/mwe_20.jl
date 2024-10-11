# generate test cases for the calculation of roll, pitch and yaw
using LinearAlgebra, Rotations

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

function is_right_handed_orthonormal(x, y, z)
    R = [x y z]
    R*R' ≈ I && det(R) ≈ 1
end

@assert is_right_handed_orthonormal(x, y, z)

# R = Yaw * Pitch * Roll

yaw = deg2rad(20)
Yaw = AngleAxis(yaw, z[1], z[2], z[3])
x1 = Yaw * x
y1 = Yaw * y
z1 = Yaw * z
println("Yaw: ", rad2deg(yaw), " x1: ", x1, " y1: ", y1, " z1: ", z1)
pitch = deg2rad(30)
Pitch = AngleAxis(pitch, y1[1], y1[2], y1[3])
x2 = Pitch * x
y2 = Pitch * y
z2 = Pitch * z
println("Pitch: ", rad2deg(pitch), " x2: ", x2, " y2: ", y2, " z2: ", z2)
R = Yaw * Pitch
x3 = R * x
y3 = R * y
z3 = R * z
println("Yaw: ", rad2deg(yaw), ", Pitch: ", rad2deg(pitch), " x3: ", x3, " y3: ", y3, " z3: ", z3)

