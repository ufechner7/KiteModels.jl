using Pkg
using Test, BenchmarkTools, StaticArrays, LinearAlgebra, KiteUtils
using KiteModels, KitePodModels

set_data_path(joinpath(dirname(dirname(pathof(KiteModels))), "data"))
set = deepcopy(load_settings("system.yaml"))
kcu::KCU = KCU(set)
kps5::KPS5 = KPS5(kcu)

pos, vel = nothing, nothing
poss, vels = nothing, nothing

@testset verbose = true "KPS5 tests...." begin

function set_defaults()
    KiteModels.clear!(kps5)
    kps5.set.l_tether = 150.0
    kps5.set.elevation = 60.0
    kps5.set.area = 20.0
    kps5.set.rel_side_area = 50.0
    kps5.set.v_wind = 8.0
    kps5.set.mass = 11.4
    kps5.set.damping =  2 * 473.0
    kps5.set.alpha = 1.0/7
    kps5.set.c_s = 0.6
    kps5.set.kcu_diameter = 0
end

function init_392()
    KiteModels.clear!(kps5)
    kps5.set.l_tether = 392.0
    kps5.set.elevation = 70.0
    kps5.set.area = 10.0
    kps5.set.rel_side_area = 50.0
    kps5.set.v_wind = 9.1
    kps5.set.mass = 6.2
    kps5.set.c_s = 0.6
end

function init_150()
    KiteModels.clear!(kps5)
    kps5.set.l_tether = 150.0
    kps5.set.elevation = 70.0
    kps5.set.area = 10.18
    kps5.set.rel_side_area = 30.6
    kps5.set.v_wind = 9.1
    kps5.set.mass = 6.21
    kps5.set.c_s = 0.6
    kps5.set.damping = 473.0     # unit damping
    kps5.set.c_spring = 614600.0 # unit spring coefficent
    kps5.set.width = 4.9622
end

function init3()
    kps5.set.alpha =  0.08163
    KiteModels.clear!(kps5)
    kps5.set.l_tether = 150.0 # - kps5.set.height_k - kps5.set.h_bridle
    kps5.set.area = 10.18
    kps5.set.rel_side_area = 30.6
    kps5.set.mass = 6.21
    kps5.set.c_s = 0.6
    kps5.set.damping = 473.0     # unit damping coefficient
    kps5.set.c_spring = 614600.0 # unit spring coefficent
    kps5.set.width = 4.9622
    kps5.set.elevation = 70.7 
    kps5.set.profile_law = Int(EXPLOG)
    pos, vel = KiteModels.init_pos_vel(kps5)
    posd = copy(vel)
    veld = zero(vel)
    height = 134.14733504839947
    kps5.v_wind .= kps5.v_wind_gnd * calc_wind_factor(kps5.am, height)
    kps5.stiffness_factor = 1.0
    KiteModels.init_springs!(kps5)
    return pos, vel, posd, veld
end

set_defaults()

@testset "calc_rho              " begin
    @test isapprox(calc_rho(kps5.am, 0.0), 1.225, atol=1e-5) 
    @test isapprox(calc_rho(kps5.am, 100.0), 1.210756, atol=1e-5) 
    kps5_= KPS5(kcu)
    @test kps5_ isa KPS5
end

end
nothing