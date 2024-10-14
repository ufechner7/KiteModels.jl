# find the reverse function of the given function
using Rotations, StaticArrays, Pkg, Test
pkg"add KiteUtils#main"
using KiteUtils

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

function frame2quat(x, y, z)
    return QuatRotation(rot3d(x, y, z))
end

quat2viewer(orient::AbstractVector) = quat2viewer(QuatRotation(orient))
function calc_orient_quat_(x, y, z)  
    x = enu2ned(x)
    y = enu2ned(y)
    z = enu2ned(z)            

    ax = @SVector [1, 0, 0]
    ay = @SVector [0, 1, 0]
    az = @SVector [0, 0, 1]
    rotation = rot3d(ax, ay, az, x, y, z)
    q = QuatRotation(rotation)
    return Rotations.params(q)
end

quat2frame_(orient::AbstractVector) = quat2frame_(QuatRotation(orient))
function quat2frame_(q::QuatRotation)
    x = [0,  1.0, 0]
    y = [1.0,  0, 0]
    z = [0,    0, -1.0]
    return q*x, q*y, q*z
end

@testset "Testing quat2frame_ ...." begin
    global q
    x = [1, 0, 0]
    y = [0, 1, 0]
    z = [0, 0, 1]
    q = calc_orient_quat_(x, y, z)
    x1, y1, z1 = quat2frame_(q)
    @test x1 ≈ x
    @test y1 ≈ y
    @test z1 ≈ z
end
nothing