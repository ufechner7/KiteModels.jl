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

yaw = deg2rad(20)
Yaw = RFR.EulerAngleAxis(yaw, z)

pitch = deg2rad(30)
Pitch = RFR.EulerAngleAxis(pitch, y)

roll = deg2rad(40)
Roll = RFR.EulerAngleAxis(roll, x)

D= RFR.angleaxis_to_dcm(Yaw)
euler = RFR.dcm_to_angle(D, :ZYX)
yaw = -euler.a1
println("yaw: ", rad2deg(yaw))

D= RFR.angleaxis_to_dcm(Pitch)
euler = RFR.dcm_to_angle(D, :ZYX)
pitch = euler.a3
println("pitch: ", rad2deg(pitch))

D= RFR.angleaxis_to_dcm(Roll)
euler = RFR.dcm_to_angle(D, :ZYX)
roll = euler.a2
println("roll: ", rad2deg(roll))

R = Pitch * Roll * Yaw
D = RFR.angleaxis_to_dcm(R)
euler = RFR.dcm_to_angle(D, :ZYX)
yaw = -euler.a1
println("yaw: ", rad2deg(yaw))
pitch = euler.a3
println("pitch: ", rad2deg(pitch))
roll = euler.a2
println("roll: ", rad2deg(roll))