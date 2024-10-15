# generate test cases for the calculation of roll, pitch and yaw
# convert all vectors to NED reference frame first
using LinearAlgebra
import ReferenceFrameRotations as RFR

# z-y′-x″ (intrinsic rotations) or x-y-z (extrinsic rotations): 
# the intrinsic rotations are known as: yaw, pitch and roll

# x: from trailing edge to leading edge
# y: to the right looking in flight direction
# z: down

yaw = deg2rad(-90)
pitch = deg2rad(0)
roll = deg2rad(0)

function enu2ned(vec::AbstractVector)  
    R = [0 1 0; 1 0 0; 0 0 -1]
    R*vec
end

function ned2enu(vec::AbstractVector)  
    R = [0 1 0; 1 0 0; 0 0 -1]
    R*vec
end

function euler2rot(roll, pitch, yaw)
    φ      = roll
    R_x = [1    0       0;
              0  cos(φ) -sin(φ);
              0  sin(φ)  cos(φ)]
    θ      = pitch          
    R_y = [ cos(θ)  0  sin(θ);
                 0     1     0;
              -sin(θ)  0  cos(θ)]
    ψ      = yaw
    R_z = [cos(ψ) -sin(ψ) 0;
              sin(ψ)  cos(ψ) 0;
                 0       0   1]
    R   = R_z * R_y * R_x
    return R
end

function rot2euler(rot)  
    D = rot
    pitch = asin(−D[3,1])
    roll  = atan(D[3,2], D[3,3])
    yaw   = atan(D[2,1], D[1,1])
    return roll, pitch, yaw
end

# If kite (x axis) is pointing to the north, and is at zenith, then in ENUs reference frame:
# - x = 0, 1, 0
# - y = 1, 0, 0
# - z = 0, 0,-1
# This would be equal to the NED reference frame.
x = enu2ned([0, 1, 0])
y = enu2ned([1, 0, 0])
z = enu2ned([0, 0,-1])

# D1 = RFR.angle_to_dcm(yaw, pitch, roll, :ZYX)
D1 = euler2rot(roll, pitch, yaw)
x4 = D1 * x
y4 = D1 * y
z4 = D1 * z
println("Yaw: ", rad2deg(yaw), ", Pitch: ", rad2deg(pitch), ", Roll: ", rad2deg(roll), 
        "\nx = ", ned2enu(x4), "\ny = ", ned2enu(y4), "\nz = ", ned2enu(z4))

roll, pitch, yaw = rot2euler(D1)
println("Yaw: ", rad2deg(yaw), ", Pitch: ", rad2deg(pitch), ", Roll: ", rad2deg(roll))
