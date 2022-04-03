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
    kps4.set.area = 10.18
    kps4.set.rel_side_area = 30.6
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
    @test sp[8].length ≈ 2.4403799411263991
    @test sp[9].length ≈  3.3353120022370204
    @test sp[10].length ≈ 3.3353120022370204
    @test sp[11].length ≈ 5.4912468596622306
    @test sp[11].c_spring ≈ 49165.631334315847
    @test sp[12].length ≈ 5.4912468596622306
    # TODO also test spring 13 .. 15
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
        kps4.v_wind_tether .= KVec3(8.0, 0.1, 0.0)
        bytes = @allocated KiteModels.calc_particle_forces(kps4, pos1, pos2, vel1, vel2, spring, stiffnes_factor, se().segments, se().d_tether/1000.0, rho, i)
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
    kps4.set.elevation = 60.0
    pos, vel = KiteModels.init_pos_vel(kps4)
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

@testset "inner_loop            " begin
    kps4.set.alpha = 1.0/7.0
    init_150()
    kps4.set.elevation = 60.0
    kps4.set.profile_law = Int(EXP)
    for i in 1:se().segments + KiteModels.KITE_PARTICLES + 1 
        kps4.forces[i] .= zeros(3)
    end
    pos, vel = KiteModels.init_pos_vel(kps4)
    v_wind_gnd = KVec3(7.0, 0.1, 0.0)
    stiffnes_factor = 0.5
    segments = kps4.set.segments
    d_tether = kps4.set.d_tether/1000.0
    KiteModels.inner_loop(kps4, pos, vel, v_wind_gnd, stiffnes_factor, segments, d_tether)
    forces =  [[ -1.1039795506035208  -0.0210281466470539   0.6374018106640786]
               [ -2.6112444243501161  -0.0497374193345597   1.5075837513184493]
               [ -3.2469329482256093  -0.0618458809786162   1.8746176116987205]
               [ -3.6500726231166869  -0.0695247045888339   2.1073704115748102]
               [ -3.9578038144544192  -0.0753862315039757   2.2850390976806585]
               [ -4.2101657193446407  -0.0801931096212583   2.4307403112207719]
               [-36.5076651709692399  -0.0567130032327103 -65.4202414843858833]
               [-31.8997735210856632  -0.0144731948993382  20.0026639693058534]
               [ 39.0565328469252862  -0.0092205388684263  47.4026845127173857]
               [  9.7854297156670409  84.4721933009509485   0.5027899265154456]
               [  9.7877962918810688 -84.5072746928884726   0.5210774922270867]]
    for i in 1:se().segments + KiteModels.KITE_PARTICLES + 1
        @test all(forces[i,:] .≈ kps4.forces[i])
    end
end

@testset "calc_aero_forces      " begin
    kps4.set.alpha = 1.0/7.0
    init_150()
    kps4.set.elevation = 60.0
    kps4.set.profile_law = Int(EXP)
    for i in 1:se().segments + KiteModels.KITE_PARTICLES + 1 
        kps4.forces[i] .= zeros(3)
    end
    pos, vel = KiteModels.init_pos_vel(kps4)
    rho = 1.25
    kps4.v_wind .= KVec3(8.0, 0.2, 0.0)
    alpha_depower = 0.1
    rel_steering = -0.1
    kps4.set.alpha_zero = 5.0
    for i in 1:se().segments + KiteModels.KITE_PARTICLES + 1 
        kps4.forces[i] .= zeros(3)
    end
    KiteModels.calc_aero_forces(kps4, pos, vel, rho, alpha_depower, rel_steering)
    forces = [[   0.                    0.                    0.                ]
              [   0.                    0.                    0.                ]
              [   0.                    0.                    0.                ]
              [   0.                    0.                    0.                ]
              [   0.                    0.                    0.                ]
              [   0.                    0.                    0.                ]
              [   0.                    0.                    0.                ]
              [   0.                    0.                    0.                ]
              [-179.1688872511660122   -4.4791993800719236 -308.8002807236504736]
              [ -11.5996366068210843  -82.7996046395915926   -1.1901722621694235]
              [ -11.3939179308052196   53.1539392663951702    0.7708381078373823]]
    for i in 1:se().segments + KiteModels.KITE_PARTICLES + 1
        @test all(forces[i,:] .≈ kps4.forces[i])
    end
end

function init2()
    kps4.set.alpha = 1.0/7.0
    init_150()
    kps4.set.elevation = 60.0
    kps4.set.profile_law = Int(EXP)
    pos, vel = KiteModels.init_pos_vel(kps4)
    posd = copy(vel)
    veld = zero(vel)
    kps4.v_wind_gnd .= [7.0, 0.1, 0.0]
    height = 134.14733504839947
    kps4.v_wind .= kps4.v_wind_gnd * calc_wind_factor(kps4, height)
    kps4.stiffness_factor = 0.5
    kps4.set.alpha = 1.0/7.0
    length = 150.0
    kps4.segment_length = length/se().segments
    for i in 1:se().segments + KiteModels.KITE_PARTICLES + 1 
        kps4.forces[i] .= zeros(3)
    end
    KiteModels.init_springs(kps4)
    return pos, vel, posd, veld
end

@testset "test_loop          " begin
    init2()
    pos, vel, posd, veld = init2()
    KiteModels.loop(kps4, pos, vel, posd, veld)
    res1=[[-0.        0.000001 -0.      ]
          [ 0.        0.        0.      ]
          [ 0.        0.        0.      ]
          [ 0.        0.        0.      ]
          [ 0.        0.        0.      ]
          [ 0.        0.        0.      ]
          [ 0.        0.        0.      ]
          [ 0.        0.        0.      ]
          [ 0.        0.        0.      ]
          [ 0.        0.        0.      ]
          [ 0.        0.        0.      ]]
    for i in 1:se().segments + KiteModels.KITE_PARTICLES + 1 
        @test all(kps4.res1[i] .≈ res1[i,:])
    end
    forces=[[ -1.1039795506226362  -0.0210281466470539   0.6374018106309699]
            [ -2.6112444243501161  -0.0497374193345597   1.5075837513184491]
            [ -3.2469329482256093  -0.0618458809786162   1.8746176116987205]
            [ -3.6500726230975715  -0.0695247045888339   2.1073704116079188]
            [ -3.9578038144525074  -0.0753862315039757   2.2850390976839696]
            [ -4.2101657194038982  -0.0801931096212583   2.4307403111181349]
            [-36.507665170911892   -0.0567130032327103 -65.4202414842865778]
            [-31.8997735210856632  -0.0144731948993382  20.0026639693058534]
            [ 39.0565328469252862  -0.0092205388684263  47.4026845127173857]
            [  9.7854297156670409  84.4721933009509485   0.5027899265154456]
            [  9.7877962918810688 -84.5072746928884726   0.5210774922270867]]
    for i in 1:se().segments + KiteModels.KITE_PARTICLES + 1 
        @test all(forces[i,:] .≈ kps4.forces[i])
    end
    res2 = [[  0.000001             0.000001             0.000001          ]
            [ -9.4954342703640595  -0.1808633430347626  15.2921227320670887]
            [-11.8070289026385815  -0.2248941126495136  16.6267913152680755]
            [-13.2729913567184443  -0.2528171075957596  17.4731651331197071]
            [-14.3920138707363918  -0.2741317509235479  18.1192330824871632]
            [-15.3096935251050859  -0.2916113077136664  18.6490556767932176]
            [ -4.2761540463732821  -0.0066428115060276   2.1473070003763892]
            [-10.9294458221419344  -0.0049587812722576  16.663278503890723 ]
            [ 29.6664941261243911  -0.0070037210740637  45.816049670887935 ]
            [  9.9103998578748431  85.5509913012598417  10.3192110782116959]
            [  9.9127966577351092 -85.5865207191570505  10.3377321952086678]]
    for i in 1:se().segments + KiteModels.KITE_PARTICLES + 1
        @test all(res2[i,:] .≈ kps4.res2[i])
    end
end

@testset "test_residual!       " begin
    init2()
    kps4.alpha_depower = -0.820659579962 
    kps4.stiffness_factor = 0.04
    kps4.set.alpha_zero = 0.0
    res =  zeros(MVector{6*(kps4.set.segments+4)+2, SimFloat})
    y0, yd0 = KiteModels.init_flat(kps4; old=true)
    y0s = """-0.             0.000001             -0.                   12.5000000000000036
            0.000001             21.6506350946109656   25.0000000000000071
            0.000001             43.3012701892219312   37.5000000000000071
            0.000001             64.9519052838329003   50.0000000000000142
            0.000001             86.6025403784438623   62.5000000000000142
            0.000001            108.2531754730548244   75.0000000000000142
            0.000001            129.9038105676658006   76.5905217485472747    0.
        134.6435550484500823   78.5650000003636109    0.                  136.0785716987731178
        77.4500000002498865    2.4810999999999996  134.1473350483994693
        77.4500000002498865   -2.4810999999999996  134.1473350483994693
            0.000001              0.000001              0.000001              0.000001
            0.000001             -0.                    0.000001              0.000001
           -0.                    0.000001              0.000001             -0.
            0.000001              0.000001             -0.                    0.000001
            0.000001             -0.                    0.000001              0.000001
           -0.                    0.000001              0.000001              0.000001
            0.000001              0.000001              0.000001              0.000001
            0.000001              0.000001              0.000001              0.000001
            0.000001            150.                    0."""
    y0s = split(replace(y0s, r"[\s]+" => ","),",")
    @test length(y0s) == length(y0)
    y0_ = Float64[]
    for i in 1:length(y0s)
        push!(y0_, parse(Float64, y0s[i]))
        if ! (y0_[i] ≈ y0[i])
            println(y0_[i]," != ", y0[i])
        end
    end
    @test sum(y0) ≈ 1716.23568958026658
    @test sum(yd0) ≈ -98.09994900000005
    time = 0.0
    @test length(res) == length(y0)-6
    residual!(res, yd0, y0, kps4, time)
    println("length(res): $(length(res))")
    res1, res2 = kps4.res1, kps4.res2
    height = calc_height(kps4)
    # @test height ≈ 134.14733504839947
    res1_= zeros(MVector{(kps4.set.segments+4+1), KVec3})
    res1_[1][2]=1e-6
    @test res1 == res1_
    res2_= [[   0.000001              0.000001              0.000001          ]
            [  -9.4954332703640603   -0.1808623430347626    5.4821227320670882]
            [ -11.8070279026385823   -0.2248931126495136    6.816791315268075 ]
            [ -13.2729903567184451   -0.2528161075957597    7.6631651331197066]
            [ -14.3920128707363926   -0.2741307509235479    8.3092330824871627]
            [ -15.3096925251050866   -0.2916103077136664    8.8390556767932171]
            [  -0.7543438292598375   -0.0077244009806506   -0.407464478647972 ]
            [  -1.2766764661730894   -0.0049560787501444    0.6267331714897555]
            [-170.1483123967530844   -2.4632143900241044  327.4738349406080147]
            [ -19.5074216613609686 -108.106358305295231    -0.6312062960364813]
            [ -22.7934890811423543  107.5043152723890785    1.2645963793717296]]
    @test res2_[1, :] == res2[1]
    @test all(res2_[2, :] .≈ res2[2])
    for i in 1:se().segments + KiteModels.KITE_PARTICLES + 1
        @test all(res2_[i, :] .≈ res2[i])
    end

end

@testset "test_find_steady_state" begin
    KiteModels.set_depower_steering(kps4, 0.25, 0.0)
    init2()
    kps4.alpha_depower = -0.820659579962 
    kps4.v_wind_gnd .= [9.0, 0.0, 0.0]
    height = 134.14733504839947
    kps4.v_wind .= kps4.v_wind_gnd * calc_wind_factor(kps4, height)    
    kps4.stiffness_factor = 1.0
    kps4.set.alpha_zero = 0.0   
    res1, res2 = find_steady_state(kps4, true) 
#    @test norm(res2) < 1e-5                            # velocity and acceleration must be near zero
#    pre_tension = KiteModels.calc_pre_tension(kps)
#    @test pre_tension > 1.0001
#    @test pre_tension < 1.01
#    @test unstretched_length(kps) ≈ 392.0              # initial, unstreched tether lenght
#    @test tether_length(kps) ≈ 392.1861381318156 # real, streched tether length
#    @test winch_force(kps) ≈ 276.25776695110034        # initial force at the winch [N]
#    lift, drag = lift_drag(kps)
#    @test lift ≈ 443.63303000106197                    # initial lift force of the kite [N]
#    @test drag ≈ 94.25223134952152                     # initial drag force of the kite [N]
#    @test lift_over_drag(kps) ≈ 4.706870316479931      # initlal lift-over-drag
#    @test norm(v_wind_kite(kps)) ≈ 9.107670173739065   # inital wind speed at the height of the kite [m/s]
end

end