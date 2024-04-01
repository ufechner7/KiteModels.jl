using Test, BenchmarkTools, StaticArrays, LinearAlgebra, KiteUtils
using KiteModels, KitePodModels

if ! @isdefined kcu
    const kcu = KCU(se())
end
if ! @isdefined kps4
    const kps4 = KPS4(kcu)
end

msg = String[]
@testset verbose = true "KPS4 benchmarking...." begin

function set_defaults()
    KiteModels.clear!(kps4)
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
    KiteModels.clear!(kps4)
    kps4.set.l_tether = 392.0
    kps4.set.elevation = 70.0
    kps4.set.area = 10.0
    kps4.set.rel_side_area = 50.0
    kps4.set.v_wind = 9.1
    kps4.set.mass = 6.2
    kps4.set.c_s = 0.6
end

function init_150()
    KiteModels.clear!(kps4)
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

function init2()
    kps4.set.alpha = 1.0/7.0
    init_150()
    kps4.set.elevation = 60.0
    kps4.set.profile_law = Int(FAST_EXP)
    pos, vel = KiteModels.init_inner(kps4)
    posd = copy(vel)
    veld = zero(vel)
    kps4.v_wind_gnd .= [7.0, 0.1, 0.0]
    kps4.stiffness_factor = 0.5
    kps4.set.alpha = 1.0/7.0
    length = 150.0
    kps4.segment_length = length/se().segments
    for i in 1:se().segments + KiteModels.KITE_PARTICLES + 1 
        kps4.forces[i] .= zeros(3)
    end
    return pos, vel, posd, veld
end

set_defaults()

# benchmark calc_particle_forces!
t = @benchmark KiteModels.calc_particle_forces!(kps4, pos1, pos2, vel1, vel2, spring, segments, d_tether, rho, i) setup=(pos1 = KVec3(1.0, 2.0, 3.0);  
                                        pos2 = KVec3(2.0, 3.0, 4.0); vel1 = KVec3(3.0, 4.0, 5.0); vel2 = KVec3(4.0, 5.0, 6.0); kps4.v_wind_tether.=KVec3(8.0, 0.1, 0.0); spring=kps4.springs[1];
                                        kps4.stiffness_factor = 0.5; segments=6.0; d_tether=se().d_tether/1000.0; rho=kps4.set.rho_0; i=rand(1:se().segments + KiteModels.KITE_PARTICLES + 1))
@test t.memory == 0
global msg
push!(msg, ("Mean time calc_particle_forces!: $(round(mean(t.times), digits=1)) ns"))

# benchmark inner_loop!
pos, vel = KiteModels.init_inner(kps4)
t = @benchmark KiteModels.inner_loop!(kps4, pos, vel, v_wind_gnd, segments, d_tether) setup=(kps4.set.elevation = 60.0; kps4.set.profile_law = Int(FAST_EXP);
                                      kps4.set.alpha = 1.0/7.0; pos = $pos; vel=$vel;
                                      v_wind_gnd = KVec3(7.0, 0.1, 0.0); kps4.stiffness_factor = 0.5; segments = kps4.set.segments; d_tether = kps4.set.d_tether/1000.0)
push!(msg, ("Mean time inner_loop!:          $(round(mean(t.times), digits=1)) ns"))
@test t.memory == 0

# benchmark calc_aero_forces!
kps4.set.alpha = 1.0/7.0
init_150()
kps4.set.elevation = 60.0
kps4.set.profile_law = Int(FAST_EXP)
for i in 1:se().segments + KiteModels.KITE_PARTICLES + 1 
    kps4.forces[i] .= zeros(3)
end
pos, vel = KiteModels.init_inner(kps4)
rho = 1.25
kps4.v_wind .= KVec3(8.0, 0.2, 0.0)
alpha_depower = 0.1
rel_steering = -0.1
kps4.set.alpha_zero = 5.0
for i in 1:se().segments + KiteModels.KITE_PARTICLES + 1 
    kps4.forces[i] .= zeros(3)
end
t = @benchmark KiteModels.calc_aero_forces!($kps4, $pos, $vel, $rho, $alpha_depower, $rel_steering)
push!(msg, ("Mean time calc_aero_forces!:    $(round(mean(t.times), digits=1)) ns"))
@test t.memory == 0

# benchmark loop!
init2()
pos, vel, posd, veld = init2()
t = @benchmark KiteModels.loop!($kps4, $pos, $vel, $posd, $veld)
push!(msg, ("Mean time loop!:                $(round(mean(t.times), digits=1)) ns"))
@test t.memory == 0

# benchmark residual!
init2()
kps4.alpha_depower = -0.820659579962 
kps4.stiffness_factor = 0.04
kps4.set.alpha_zero = 0.0
res =  zeros(MVector{6*(kps4.set.segments+KiteModels.KITE_PARTICLES), SimFloat})
y0, yd0 = KiteModels.init(kps4)
time = 0.0
t = @benchmark residual!($res, $yd0, $y0, $kps4, $time)
push!(msg, ("Mean time residual!:           $(round(mean(t.times), digits=1)) ns"))
@test t.memory == 0
# time using Python/ Numba: 8.94 µs, time using Julia 1.7.2: 1.6µs, Julia 1.8.0: 1.244µs
# Julia 1.9 on Ryzen:  816.1 ns
# Julia 1.10 on Ryzen: 787.0 ns 6000 RAM
end
for i in eachindex(msg)
    println(msg[i])
end

