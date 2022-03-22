using Test, BenchmarkTools, StaticArrays, LinearAlgebra, KiteUtils
using KiteModels, KitePodModels

if ! @isdefined kcu
    const kcu = KCU()
end
if ! @isdefined kps4
    const kps4 = KPS4(kcu)
end

msg = String[]
@testset verbose = true "KPS4 benchmarking...." begin

function set_defaults()
    KiteModels.clear(kps4)
    kps4.set.l_tether = 150.0
    kps4.set.elevation = 60.0
    kps4.set.area = 20.0
    kps4.set.rel_side_area = 50.0
    kps4.set.v_wind = 8.0
    kps4.set.mass = 11.4
    kps4.set.damping =  2 * 473.0
    kps4.set.alpha = 1.0/7
    kps4.set.c_s = 0.6
    
end

function init_392()
    KiteModels.clear(kps4)
    kps4.set.l_tether = 392.0
    kps4.set.elevation = 70.0
    kps4.set.area = 10.0
    kps4.set.rel_side_area = 50.0
    kps4.set.v_wind = 9.1
    kps4.set.mass = 6.2
    kps4.set.c_s = 0.6
end

function init_150()
    KiteModels.clear(kps4)
    kps4.set.l_tether = 150.0
    kps4.set.elevation = 70.0
    kps4.set.area = 10.0
    kps4.set.rel_side_area = 50.0
    kps4.set.v_wind = 9.1
    kps4.set.mass = 6.21
    kps4.set.c_s = 0.6
    kps4.set.damping = 473.0     # unit damping coefficient
    kps4.set.c_spring = 614600.0 # unit spring coefficent
    kps4.set.width = 4.9622
end

set_defaults()

@testset "calc_particle_forces  " begin
    init_150()
    sp = KiteModels.init_springs(kps4)
    pos1 = [1.0, 2.0, 3.0]
    pos2 = [2.0, 3.0, 4.0]
    vel1 = [3.0, 4.0, 5.0]
    vel2 = [4.0, 5.0, 6.0]
    rho = kps4.set.rho_0
    for i in 1:se().segments + KiteModels.KITE_PARTICLES + 1 
        kps4.forces[i] .= zeros(3)
    end
    for i in 1:length(kps4.springs)
        spring = kps4.springs[i]
        stiffnes_factor = 0.5
        kps4.v_wind_tether .= [8.0, 0.1, 0.0]
        KiteModels.calc_particle_forces(kps4, pos1, pos2, vel1, vel2, spring, stiffnes_factor, se().segments, se().d_tether/1000.0, rho, i)
    end
    # Python output
    res=[[ 18550.4729309395152086  18550.6132232745367219  18550.6305627766196267]
         [    -0.1986161147506209      0.0819685552924057      0.1166475594582153]
         [    -0.1986161147506209      0.0819685552924057      0.1166475594582153]
         [    -0.1986161147506209      0.0819685552924057      0.1166475594582153]
         [    -0.1986161147506209      0.0819685552924057      0.1166475594582153]
         [    -0.1986161147506209      0.0819685552924057      0.1166475594582153]
         [-32417.6381463687685027 -32417.0769770286824496 -32417.0076190203544684]
         [-20528.0512582440096594 -20527.4900889039272442 -20527.420730895599263 ]
         [ 12986.35257788861054    12986.7734548936750798  12986.8254733999201562]
         [ 23289.9810739697131794  23290.5422433097955945  23290.6116013181235758]
         [ -1883.1033393325606085  -1882.5421699924754648  -1882.4728119841502121]]
    for i in 1:se().segments + KiteModels.KITE_PARTICLES + 1
        @test all(res[i,:] .≈ kps4.forces[i])
    end
end

t = @benchmark KiteModels.calc_particle_forces(kps4, pos1, pos2, vel1, vel2, spring, stiffnes_factor, segments, d_tether, rho, i) setup=(pos1 = KVec3(1.0, 2.0, 3.0);  
                                        pos2 = KVec3(2.0, 3.0, 4.0); vel1 = KVec3(3.0, 4.0, 5.0); vel2 = KVec3(4.0, 5.0, 6.0); kps4.v_wind_tether.=KVec3(8.0, 0.1, 0.0); spring=kps4.springs[1];
                                        stiffnes_factor = 0.5; segments=6.0; d_tether=se().d_tether/1000.0; rho=kps4.set.rho_0; i=rand(1:se().segments + KiteModels.KITE_PARTICLES + 1))
@test t.memory == 0
global msg
push!(msg, ("Mean time calc_particle_forces: $(round(mean(t.times), digits=1)) ns"))

pos, vel = KiteModels.init(kps4)
t = @benchmark KiteModels.inner_loop2(kps4, pos, vel, v_wind_gnd, stiffnes_factor, segments, d_tether) setup=(kps4.set.elevation = 60.0; kps4.set.profile_law = 1;
                                      kps4.set.alpha = 1.0/7.0; pos = $pos; vel=$vel;
                                      v_wind_gnd = KVec3(7.0, 0.1, 0.0); stiffnes_factor = 0.5; segments = kps4.set.segments; d_tether = kps4.set.d_tether/1000.0)
push!(msg, ("Mean time inner_loop2: $(round(mean(t.times), digits=1)) ns"))
@test t.memory == 0

end
println(msg[1])
println(msg[2])

