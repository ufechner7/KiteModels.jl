# find the reverse function of the given function
using Rotations, StaticArrays, Pkg, Test, LinearAlgebra
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
    @assert is_right_handed_orthonormal(x, y, z)
    R = [x y z]
    return (QuatRotation(R))
end
pure_quat2frame(q::AbstractVector) = pure_quat2frame(QuatRotation(q))
function pure_quat2frame(q::QuatRotation)
    x = enu2ned(q[1,:]) .* [1, 1, -1]
    y = enu2ned(q[2,:])
    z = enu2ned(q[3,:]) .* [1, -1, 1]
    return x, y, z
end

calc_orient_quat_(orient::AbstractVector) = calc_orient_quat_(QuatRotation(orient))
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
    global x1, y1, z1, q, rotation
    x = [1, 0, 0]
    y = [0, 1, 0]
    z = [0, 0, 1]
    q = calc_orient_quat_(x, y, z)
    x1, y1, z1 = quat2frame_(q)
    @test x1 ≈ x
    @test y1 ≈ y
    @test z1 ≈ z
    # test with a different orientation
    rotation = euler2rot(1, 0, 0)
    q1 = QuatRotation(rotation)
    q = frame2quat(x, y, z)
    x2, y2, z2 = quat2frame_(q*q1)
    q = calc_orient_quat_(x2, y2, z2)
    x3, y3, z3 = pure_quat2frame(q)
    @test x3 ≈ x2
    @test y3 ≈ y2
    @test z3 ≈ z2
end
function test(roll=0)
    rotation = euler2rot(roll, 0, 0)
    ax = @SVector [1, 0, 0]
    ay = @SVector [0, 1, 0]
    az = @SVector [0, 0, 1]
    @assert is_right_handed_orthonormal(ax, ay, az)
    q = frame2quat(ax, ay, az)
    q1 = QuatRotation(rotation)
    return q*q1
end
nothing
