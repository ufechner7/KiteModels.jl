using Test, BenchmarkTools, StaticArrays, LinearAlgebra, KiteUtils
using KiteModels, KitePodModels

if ! @isdefined kcu
    const kcu = KCU()
end
if ! @isdefined kps4
    const kps4 = KPS4(kcu)
end

@testset verbose = true "KPS4 tests...." begin

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

@testset "calc_rho              " begin
    @test isapprox(calc_rho(kps4, 0.0), 1.225, atol=1e-5) 
    @test isapprox(calc_rho(kps4, 100.0), 1.210756, atol=1e-5) 
end

@testset "initial_kite_ref_frame" begin
    vec_c    = [-15., 0., -25.98076211]
    v_app    = [10.4855, 0, -3.08324]
    x, y, z = KiteModels.initial_kite_ref_frame(vec_c, v_app)
    @test x == [-0.8660254037549957, 0.0, 0.5000000000509968]
    @test y == [0.0, 1.0, 0.0]
    @test z == [-0.5000000000509968, 0.0, -0.8660254037549957]
end

@testset "get_particles         " begin
    init_150()
    particles = KiteModels.get_particles(kps4.set.height_k, kps4.set.h_bridle, kps4.set.width, kps4.set.m_k)
    @test particles[1] == zeros(3)
    @test particles[2] == [  75.,              0.    ,       129.90381057]
    @test particles[3] == [76.590521748547275, 0.    , 134.64355504845008]
    @test particles[4] == [78.565000000363611, 0.    , 136.07857169877312]
    @test particles[5] == [77.450000000249887, 2.4811, 134.14733504839947]
    @test particles[6] == [77.450000000249887,-2.4811, 134.14733504839947]
end

@testset "init_springs          " begin
    init_150()
    sp = KiteModels.init_springs(kps4)
    # test springs
    @test length(sp) == 6 + KiteModels.KITE_SPRINGS
    @test sp[1].p1 == 1
    @test sp[1].p2 == 2
    for i in 1:6
        @test sp[i].length   ≈ 150.0/6
        @test sp[i].c_spring ≈ 2.76460154e+04
        @test sp[i].damping  ≈ 1.89200000e+01
    end
    @test sp[6].p1 == 6
    @test sp[6].p2 == 7
    @test sp[7].p1 == 7
    @test sp[7].p2 == 8
    @test sp[7].length ≈ 4.998493790987047
    @test sp[7].c_spring ≈ 54012.39452466341
    @test sp[7].damping ≈ 94.62850606174248
end

@testset "init_masses           " begin
    init_150()
    m = KiteModels.init_masses(kps4)
    @test m[1] ≈ 0.1137256540599505
    for i in 2:6
        @test m[i] ≈ 0.227451308119901
    end
    @test m[7] ≈ 8.5137256540599502
    @test m[8] ≈ 2.9187
    @test m[9] ≈ 1.31652
    @test m[11] ≈ 0.98739
    @test m[11] ≈ 0.98739
end

@testset "calc_particle_forces  " begin
    init_150()
    pos1 = KVec3(1.0, 2.0, 3.0)
    pos2 = KVec3(2.0, 3.0, 4.0)
    vel1 = KVec3(3.0, 4.0, 5.0)
    vel2 = KVec3(4.0, 5.0, 6.0)
    rho = kps4.set.rho_0
    for i in 1:se().segments + KiteModels.KITE_PARTICLES + 1 
        kps4.forces[i] .= zeros(3)
    end
    bytes = 0
    for i in 1:length(kps4.springs)
        spring = kps4.springs[i]
        stiffnes_factor = 0.5
        v_wind_tether = KVec3(8.0, 0.1, 0.0)
        bytes = @allocated KiteModels.calc_particle_forces(kps4, pos1, pos2, vel1, vel2, v_wind_tether, spring, stiffnes_factor, se().segments, se().d_tether/1000.0, rho, i)
    end
    # @test bytes == 0
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

@testset "init                  " begin
    init_150()
    kps4.set.elevation=60.0
    pos, vel = KiteModels.init(kps4)
    pos1 = [[  -0.                    0.000001             -0.                ]
            [  12.5000000000000036    0.000001             21.6506350946109656]
            [  25.0000000000000071    0.000001             43.3012701892219312]
            [  37.5000000000000071    0.000001             64.9519052838329003]
            [  50.0000000000000142    0.000001             86.6025403784438623]
            [  62.5000000000000142    0.000001            108.2531754730548244]
            [  75.0000000000000142    0.000001            129.9038105676658006]
            [  76.5905217485472747    0.                  134.6435550484500823]
            [  78.5650000003636109    0.                  136.0785716987731178]
            [  77.4500000002498865    2.4810999999999996  134.1473350483994693]
            [  77.4500000002498865   -2.4810999999999996  134.1473350483994693]]
    for i in 1:se().segments + KiteModels.KITE_PARTICLES + 1
        @test all(pos1[i,:] .≈ pos[i])
    end
end

end