# generate test cases for the calculation of roll, pitch and yaw
using LinearAlgebra
import ReferenceFrameRotations as RFR

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

# R = Yaw * Pitch * Roll

yaw = deg2rad(0)
pitch = deg2rad(0)
roll = deg2rad(-90)

D1 = RFR.angle_to_dcm(yaw, pitch, roll, :ZYX)
x4 = D1 * x
y4 = D1 * y
z4 = D1 * z
println("Yaw: ", rad2deg(yaw), ", Pitch: ", rad2deg(pitch), ", Roll: ", rad2deg(roll), 
        "\nx = ", x4, "\ny = ", y4, "\nz = ", z4)

euler = RFR.dcm_to_angle(D1, :ZYX)
yaw = euler.a1
pitch = euler.a2
roll = euler.a3
println("Yaw: ", rad2deg(yaw), ", Pitch: ", rad2deg(pitch), ", Roll: ", rad2deg(roll))
