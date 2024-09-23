using KiteModels, Rotations, LinearAlgebra, StaticArrays

# kite reference frame; postions in ENU coordinates
# wind from west, nose to west, kite at zenith

z = [0, 0, -1]
y = [0, 1, 0]
x = y × z

function orient_euler(x, y, z)
    roll = atan(y[3], z[3]) - π/2
    if roll < -π/2
       roll += 2π
    end
    pitch = asin(-x[3])
    yaw = -atan(x[2], x[1]) - π/2
    if yaw < -π/2
        yaw += 2π
    end
    SVector(roll, pitch, yaw)
end

#  returns w, i, j, k
function calc_orient_quat(x, y, z)
    pos_kite_ = [0, 0, 0]
    pos_before = pos_kite_ .+ z
   
    rotation = rot(pos_kite_, pos_before, -x)
    q = QuatRotation(rotation)
    return Rotations.params(q)
end

orient = orient_euler(x, y, z)
roll, pitch, yaw = rad2deg.(orient)
println("roll: $roll, pitch: $pitch, yaw: $yaw")

q = QuatRotation(calc_orient_quat(x, y, z))
w, i, j, k = Rotations.params(q)
println("w: $w, i: $i, j: $j, k: $k")
