# TODO: add second test case
# TODO: try to implement calc_orient_rot(x, y, z; old=true) by
#       - calculating x, y, z from the quaternion
#       - calling calc_orient_rot(x, y, z; old=true)

using Rotations, KiteUtils, StaticArrays
import ReferenceFrameRotations as RFR

function calc_orient_rot(x, y, z; old=false)
    if old
        pos_kite_ = ones(3)
        pos_before = pos_kite_ .+ z
    
        rotation = rot(pos_kite_, pos_before, -x)
    else
        # reference frame for the orientation: NED (north, east, down)
        ax = [0, 1, 0] # in ENU reference frame this is pointing to the south
        ay = [1, 0, 0] # in ENU reference frame this is pointing to the west
        az = [0, 0, -1] # in ENU reference frame this is pointing down
        rotation = rot3d(ax, ay, az, x, y, z)
    end
    return rotation
end

# swap rows i and j of a, in-place
function swaprows!(A, i, j)
    for k in axes(A, 2)
        (A[i, k], A[j, k]) = (A[j, k], A[i, k])
    end
end

function new2old(q::QuatRotation)
    x, y, z = quat2frame(q)
    rot = calc_orient_rot(x, y, z; old=true)
    return QuatRotation(rot)
end
quat2frame(q::AbstractMatrix) = quat2frame(QuatRotation(q))
function quat2frame(q::QuatRotation)
    x = [0,  1.0, 0]
    y = [1.0,  0, 0]
    z = [0,    0, -1.0]
    return q*x, q*y, q*z
end
function euler2quat(;roll=0, pitch=0, yaw=0, deg=true)
    if deg
        roll = deg2rad(roll)
        pitch = deg2rad(pitch)
        yaw = deg2rad(yaw)
    end
    return QuatRotation(RFR.angle_to_dcm(yaw, pitch, roll, :ZYX))
end 

x = [0, 1.0, 0]
y = [1, 0, 0]
z = [0, 0, -1]

N = calc_orient_rot(x, y, z)
O = calc_orient_rot(x, y, z; old=true)


# q2 is the correct result
q_old= Float32[-0.4354337 0.88402545 -0.16998994; -0.06545633 0.15724015 0.98538876; 0.897838 0.4401984 -0.010602593]
q2= [-0.435434583948209 0.16998619895389078 -0.8840256869950432; -0.8978377839922781 -0.010605054960933391 0.44019864430176375; 0.06545255333205324 0.9853893773399851 0.15723783987269374] 
x1, y1, z1 = quat2frame(q_old)
x2, y2, z2 = quat2frame(q2)
N2 = calc_orient_rot(x2, y2, z2)
O2 = calc_orient_rot(x2, y2, z2; old=true)
O22 = new2old(QuatRotation(N2))
# O2 should be equal to O22, but it isn't

