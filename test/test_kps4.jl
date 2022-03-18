using Test, BenchmarkTools, StaticArrays, LinearAlgebra, KiteUtils

using KiteModels, KitePodModels

const SEGMENTS = se().segments
if ! @isdefined kcu
    const kcu = KCU()
    const kps = KPS4(kcu)
end

function set_defaults()
    KiteModels.clear(kps)
    kps.set.l_tether = 150.0
    kps.set.elevation = 60.0
    kps.set.area = 20.0
    kps.set.rel_side_area = 50.0
    kps.set.v_wind = 8.0
    kps.set.mass = 11.4
    kps.set.damping =  2 * 473.0
    kps.set.alpha = 1.0/7
    kps.set.c_s = 0.6
    
end

function init_392()
    KiteModels.clear(kps)
    kps.set.l_tether = 392.0
    kps.set.elevation = 70.0
    kps.set.area = 10.0
    kps.set.rel_side_area = 50.0
    kps.set.v_wind = 9.1
    kps.set.mass = 6.2
    kps.set.c_s = 0.6
end

set_defaults()

@testset "calc_rho             " begin
    @test isapprox(calc_rho(kps, 0.0), 1.225, atol=1e-5) 
    @test isapprox(calc_rho(kps, 100.0), 1.210756, atol=1e-5) 
end
