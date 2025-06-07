# SPDX-FileCopyrightText: 2025 Uwe Fechner
#
# SPDX-License-Identifier: MIT

# generate test cases for the calculation of roll, pitch and yaw
using LinearAlgebra, Rotations
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
Yaw = AngleAxis(yaw, z[1], z[2], z[3])

pitch = deg2rad(30)
Pitch = AngleAxis(pitch, y[1], y[2], y[3])

roll = deg2rad(40)
Roll = AngleAxis(roll, x[1], x[2], x[3])

D= RFR.DCM(Yaw)
euler = RFR.dcm_to_angle(D, :ZYX)
yaw = euler.a1
println("yaw: ", rad2deg(yaw))

D= RFR.DCM(Pitch)
euler = RFR.dcm_to_angle(D, :ZYX)
pitch = -euler.a3
println("pitch: ", rad2deg(pitch))

D= RFR.DCM(Roll)
euler = RFR.dcm_to_angle(D, :ZYX)
roll = -euler.a2
println("roll: ", rad2deg(roll))