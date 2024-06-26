using Test, BenchmarkTools, StaticArrays, LinearAlgebra, KiteUtils, ProfileView
using KiteModels, KitePodModels

if ! @isdefined kcu
    const kcu = KCU(se())
end
if ! @isdefined kps4_3l
    const kps4_3l = KPS4_3L(kcu)
end

msg = String[]
@testset verbose = true "KPS4_3L benchmarking...." begin

function set_defaults()
    KiteModels.clear!(kps4_3l)
    kps4_3l.set.segments = 6
    kps4_3l.set.l_tether = 50.0
    kps4_3l.set.elevation = 70.8
    kps4_3l.set.radius = 2.0
    kps4_3l.set.bridle_center_distance = 4.0
    kps4_3l.set.middle_length = 1.2
    kps4_3l.set.tip_length = 0.6
    kps4_3l.set.min_steering_line_distance = 1.0
    kps4_3l.set.width = 3.0
    kps4_3l.set.aero_surfaces = 10
    KiteModels.clear!(kps4_3l)
end

function init_392()
    KiteModels.clear!(kps4_3l)
    kps4_3l.set.l_tether = 392.0
    kps4_3l.set.elevation = 70.0
    kps4_3l.set.area = 10.0
    kps4_3l.set.rel_side_area = 50.0
    kps4_3l.set.v_wind = 9.1
    kps4_3l.set.mass = 6.2
    kps4_3l.set.c_s = 0.6
end

function init_50()
    KiteModels.clear!(kps4_3l)
    kps4_3l.set.l_tether = 50.0
    kps4_3l.set.elevation = 70.0
    kps4_3l.set.v_wind = 15.51
    kps4_3l.set.mass = 0.5
    kps4_3l.set.c_s = 0.6
    kps4_3l.set.damping = 473.0     # unit damping coefficient
    kps4_3l.set.c_spring = 614600.0 # unit spring coefficent
    kps4_3l.set.width = 3.0
end

function init2()
    kps4_3l.set.alpha = 1.0/7.0
    init_50()
    kps4_3l.set.elevation = 60.0
    kps4_3l.set.profile_law = Int(FAST_EXP)
    pos, vel = KiteModels.init_pos_vel(kps4_3l)
    posd = copy(vel)
    veld = zero(vel)
    kps4_3l.v_wind_gnd .= [15.51, 0.0, 0.0]
    kps4_3l.stiffness_factor = 0.5
    lengths = [50.0, norm(pos[kps4_3l.num_C]-pos[kps4_3l.num_E-2]), norm(pos[kps4_3l.num_D]-pos[kps4_3l.num_E-1])]
    kps4_3l.segment_lengths = lengths./se().segments
    for i in 1:kps4_3l.num_A
        kps4_3l.forces[i] .= zeros(3)
    end
    return pos, vel, posd, veld
end

set_defaults()

# # benchmark calc_particle_forces!
# t = @benchmark KiteModels.calc_particle_forces!(kps4_3l, pos1, pos2, vel1, vel2, spring, d_tether, rho, i) setup=(pos1 = KVec3(1.0, 2.0, 3.0);  
#                                         pos2 = KVec3(2.0, 3.0, 4.0); vel1 = KVec3(3.0, 4.0, 5.0); vel2 = KVec3(4.0, 5.0, 6.0); kps4_3l.v_wind_tether.=KVec3(8.0, 0.1, 0.0); spring=kps4_3l.springs[1];
#                                         kps4_3l.stiffness_factor = 0.5; segments=6.0; d_tether=se().d_tether/1000.0; rho=kps4_3l.set.rho_0; i=rand(1:se().segments + KiteModels.KITE_PARTICLES + 1))
# @test t.memory == 0
# global msg
# push!(msg, ("Mean time calc_particle_forces!: $(round(mean(t.times), digits=1)) ns"))

# # benchmark inner_loop!
# pos, vel = KiteModels.init_pos_vel(kps4_3l)
# t = @benchmark KiteModels.inner_loop!(kps4_3l, pos, vel, v_wind_gnd, d_tether) setup=(kps4_3l.set.elevation = 70.0; kps4_3l.set.profile_law = Int(FAST_EXP);
#                                       kps4_3l.set.alpha = 1.0/7.0; pos = $pos; vel=$vel;
#                                       v_wind_gnd = KVec3(7.0, 0.1, 0.0); kps4_3l.stiffness_factor = 0.5; segments = kps4_3l.set.segments; d_tether = kps4_3l.set.d_tether/1000.0)
# push!(msg, ("Mean time inner_loop!:          $(round(mean(t.times), digits=1)) ns"))
# @test t.memory == 0

# benchmark calc_aero_forces!
init_50()
kps4_3l.set.elevation = 70.0
kps4_3l.set.profile_law = Int(FAST_EXP)
for i in 1:se().segments + KiteModels.KITE_PARTICLES + 1 
    kps4_3l.forces[i] .= zeros(3)
end
pos, vel = KiteModels.init_inner(kps4_3l)
rho = 1.25
kps4_3l.v_wind .= KVec3(15.51, 0.0, 0.0)
for i in 1:kps4_3l.num_A
    kps4_3l.forces[i] .= zeros(3)
end
@profview [KiteModels.calc_aero_forces!(kps4_3l, pos, vel) for _ in 1:1000]
@profview [KiteModels.calc_aero_forces!(kps4_3l, pos, vel) for _ in 1:1000]
t = @benchmark KiteModels.calc_aero_forces!($kps4_3l, $pos, $vel)
println("Mean time calc_aero_forces!:    $(round(mean(t.times), digits=1)) ns")
push!(msg, ("Mean time calc_aero_forces!:    $(round(mean(t.times), digits=1)) ns"))
@test t.memory == 0
# best: 9360 == 0

# # benchmark loop!
# init2()
# pos, vel, posd, veld = init2()
# t = @benchmark KiteModels.loop!($kps4_3l, $pos, $vel, $posd, $veld)
# push!(msg, ("Mean time loop!:                $(round(mean(t.times), digits=1)) ns"))
# @test t.memory == 0

# # benchmark residual!
# init2()
# kps4_3l.stiffness_factor = 0.04
# res = zeros(MVector{6*(kps4_3l.num_A-5)+4+6, SimFloat})
# y0, yd0 = KiteModels.init(kps4_3l)
# time = 0.0
# t = @benchmark residual!($res, $yd0, $y0, $kps4_3l, $time)
# push!(msg, ("Mean time residual!:           $(round(mean(t.times), digits=1)) ns"))
# @test t.memory == 0


end
for i in eachindex(msg)
    println(msg[i])
end

