using KiteModels, Rotations, LinearAlgebra, StaticArrays

# kite reference frame; postions in ENU coordinates
# wind from west, nose to west, kite at zenith

z = [0, 0, -1]
y = [0, -1, 0]
x = y × z

function calc_orient_quat(x, y, z)
    # reference: NED
    ax = [0, 1, 0]
    ay = [1, 0, 0]
    az = [0, 0, -1]
    rotation = rot3d(ax, ay, az, x, y, z)
    q = QuatRotation(rotation)
    return Rotations.params(q)
end

function quat2euler(q)
    # Convert quaternion to RotXYZ
    rot = RotXYZ(q)
    
    # Extract roll, pitch, and yaw from RotXYZ
    roll = rot.theta1
    pitch = rot.theta2
    yaw = rot.theta3
    
    return roll, pitch, yaw
end

q = QuatRotation(calc_orient_quat(x, y, z))
w, i, j, k = Rotations.params(q)
println("w: $w, i: $i, j: $j, k: $k")

roll, pitch, yaw = rad2deg.(Rotations.params(RotXYZ(q)))
println("roll: $roll, pitch: $pitch, yaw: $yaw")

roll, pitch, yaw = rad2deg.(quat2euler(q))
println("roll: $roll, pitch: $pitch, yaw: $yaw")
