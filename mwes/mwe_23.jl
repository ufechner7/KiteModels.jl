using Rotations, KiteUtils

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

# swap columns i and j of a, in-place
function swapcols!(a::AbstractMatrix, i, j)
    for k in axes(a,1)
        a[k,i],a[k,j] = a[k,j],a[k,i]
    end
end

# swap rows i and j of a, in-place
function swaprows!(A, i, j)
    for k in axes(A, 2)
        (A[i, k], A[j, k]) = (A[j, k], A[i, k])
    end
end

function new2old(rot)
    x = rot
    swaprows!(x, 2, 3)
    x[1, :] .*= -1
    return x
end

x = [0, 1.0, 0]
y = [1, 0, 0]
z = [0, 0, -1]

N = calc_orient_rot(x, y, z)
O = calc_orient_rot(x, y, z; old=true)

