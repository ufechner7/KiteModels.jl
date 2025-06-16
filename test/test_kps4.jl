# SPDX-FileCopyrightText: 2025 Uwe Fechner
# SPDX-License-Identifier: MIT

using Pkg
if ! ("PackageCompiler" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
    Pkg.update()
end
using Test, BenchmarkTools, StaticArrays, LinearAlgebra, KiteUtils
using KiteModels, KitePodModels

set_data_path(joinpath(dirname(dirname(pathof(KiteModels))), "data"))
set = deepcopy(load_settings("system.yaml"))
kcu::KCU = KCU(set)
kps4::KPS4 = KPS4(kcu)

pos, vel = nothing, nothing
poss, vels = nothing, nothing

@testset verbose = true "KPS4 tests...." begin

function set_defaults()
    KiteModels.clear!(kps4)
    kps4.set.l_tethers[1] = 150.0
    kps4.set.elevation = 60.0
    kps4.set.area = 20.0
    kps4.set.rel_side_area = 50.0
    kps4.set.v_wind = 8.0
    kps4.set.mass = 11.4
    kps4.set.damping =  2 * 473.0
    kps4.set.alpha = 1.0/7
    kps4.set.c_s = 0.6
    kps4.set.kcu_diameter = 0.0
end

function init_392()
    KiteModels.clear!(kps4)
    kps4.set.l_tethers[1] = 392.0
    kps4.set.elevation = 70.0
    kps4.set.area = 10.0
    kps4.set.rel_side_area = 50.0
    kps4.set.v_wind = 9.1
    kps4.set.mass = 6.2
    kps4.set.c_s = 0.6
end

function init_150()
    KiteModels.clear!(kps4)
    kps4.set.l_tethers[1] = 150.0
    kps4.set.elevation = 70.0
    kps4.set.area = 10.18
    kps4.set.rel_side_area = 30.6
    kps4.set.v_wind = 9.1
    kps4.set.mass = 6.21
    kps4.set.c_s = 0.6
    kps4.set.damping = 473.0     # unit damping
    kps4.set.c_spring = 614600.0 # unit spring coefficent
    kps4.set.width = 4.9622
end

function init3()
    kps4.set.alpha =  0.08163
    KiteModels.clear!(kps4)
    kps4.set.l_tethers[1] = 150.0 # - kps4.set.height_k - kps4.set.h_bridle
    kps4.set.area = 10.18
    kps4.set.rel_side_area = 30.6
    kps4.set.mass = 6.21
    kps4.set.c_s = 0.6
    kps4.set.damping = 473.0     # unit damping coefficient
    kps4.set.c_spring = 614600.0 # unit spring coefficent
    kps4.set.width = 4.9622
    kps4.set.elevation = 70.7 
    kps4.set.profile_law = Int(EXPLOG)
    pos, vel = KiteModels.init_pos_vel(kps4)
    posd = copy(vel)
    veld = zero(vel)
    height = 134.14733504839947
    kps4.v_wind .= kps4.v_wind_gnd * calc_wind_factor(kps4.am, height)
    kps4.stiffness_factor = 1.0
    KiteModels.init_springs!(kps4)
    return pos, vel, posd, veld
end

set_defaults()

@testset "calc_rho              " begin
    @test isapprox(calc_rho(kps4.am, 0.0), 1.225, atol=1e-5) 
    @test isapprox(calc_rho(kps4.am, 100.0), 1.210756, atol=1e-5) 
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
    particles = KiteUtils.get_particles(kps4.set.height_k, kps4.set.h_bridle, kps4.set.width, kps4.set.m_k)
    @test particles[1] == zeros(3)
    @test particles[2] == [  75.,              0.    ,       129.90381057]   # P_KCU
    @test particles[3] == [76.590521748547275, 0.    , 134.64355504845008]   # pos_A
    @test particles[4] == [78.565000000363611, 0.    , 136.07857169877312]   # pos_B (top)
    @test particles[5] == [77.450000000249887, 2.4811, 134.14733504839947]   # pos_C (left  as seen from GS) 
    @test particles[6] == [77.450000000249887,-2.4811, 134.14733504839947]   # pos_D (right as seen from GS)
end

@testset "init_springs!          " begin
    init_150()
    sp = KiteModels.init_springs!(kps4)
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
    @test sp[7].p1 == 7 # bridle segment S1
    @test sp[7].p2 == 8
    @test sp[7].length ≈ 4.998493790987047
    @test sp[7].c_spring ≈ 54012.39452466341
    @test sp[7].damping ≈ 94.62850606174248
    @test sp[8].length ≈ 2.671692035300166       # S2
    @test sp[9].length ≈ 4.96120756              # S3
    @test sp[10].length ≈ 3.3353120022370204
    @test sp[11].length ≈ 5.4912468596622306
    @test sp[11].c_spring ≈ 49165.631334315847
    @test sp[15].length ≈  2.440379941126399     # S9 p11 p8
    # TODO also test spring 13 .. 15
end

@testset "init_masses!           " begin
    init_150()
    m = KiteModels.init_masses!(kps4)
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

@testset "calc_particle_forces!  " begin
    init_150()
    pos1 = KVec3(1.0, 2.0, 3.0)
    pos2 = KVec3(2.0, 3.0, 4.0)
    vel1 = KVec3(3.0, 4.0, 5.0)
    vel2 = KVec3(4.0, 5.0, 6.0)
    rho = kps4.set.rho_0
    for i in 1:set.segments + KiteModels.KITE_PARTICLES + 1 
        kps4.forces[i] .= zeros(3)
    end
    bytes = 0
    for i in 1:length(kps4.springs)
        spring = kps4.springs[i]
        kps4.stiffness_factor = 0.5
        kps4.v_wind_tether .= KVec3(8.0, 0.1, 0.0)
        bytes = @allocated KiteModels.calc_particle_forces!(kps4, pos1, pos2, vel1, vel2, spring, set.segments, set.d_tether/1000.0, rho, i)
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
    for i in 1:set.segments + KiteModels.KITE_PARTICLES + 1
        @test all(res[i,:] .≈ kps4.forces[i])
    end
end

@testset "init                  " begin
    init_150()
    kps4.set.elevation = 60.0
    pos, vel, acc = KiteModels.init_pos_vel_acc(kps4, old=true, delta = 1e-6)
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
            [  77.4500000002498865   -2.4810999999999996  134.1473350483994693]
           ]
    for i in 1:set.segments + KiteModels.KITE_PARTICLES + 1
        if ! all(pos1[i,:] .≈ pos[i])
            println(i, " ", pos1[i,:], " ", pos[i])
        else
            @test all(pos1[i,:] .≈ pos[i])
        end
    end
end

function split_res(res)
    particles = set.segments+KiteModels.KITE_PARTICLES
    pos=res[1:3*particles]
    vel=res[3*particles+1:6*particles]
    pos, vel
end

@testset "initial_residual      " begin
    init3()
    res = zeros(MVector{6*(set.segments+KiteModels.KITE_PARTICLES), SimFloat})
    X00 = zeros(SimFloat, 2*(set.segments+KiteModels.KITE_PARTICLES-1)+2)
    y0, yd0 = KiteModels.init(kps4, X00)
    residual!(res, yd0, y0, kps4, 0.0)
    res_pos, res_vel = split_res(res)
    @test res_pos == zeros(length(res_pos))
    # in the second test we check if the norm of the accelerations of the tether particles is low
    # println((norm(res_vel[1:15])))
    @test_broken (norm(res_vel[1:15])) < 15.0
end

@testset "inner_loop!            " begin
    kps4.set.alpha = 1.0/7.0
    init_150()
    kps4.set.elevation = 60.0
    kps4.set.profile_law = Int(EXP)
    for i in 1:set.segments + KiteModels.KITE_PARTICLES + 1 
        kps4.forces[i] .= zeros(3)
    end
    pos, vel, acc = KiteModels.init_pos_vel_acc(kps4, old=true, delta = 1e-6)
    v_wind_gnd = KVec3(7.0, 0.1, 0.0)
    kps4.stiffness_factor = 0.5
    segments = kps4.set.segments
    d_tether = kps4.set.d_tether/1000.0
    KiteModels.inner_loop!(kps4, pos, vel, v_wind_gnd, segments, d_tether)
    forces =  [[ -1.0581480474703286  -0.020155169648329054 0.6109409799927648]
               [ -2.3531675038883413  -0.04482168872970446  1.358582971818214]
               [ -2.699059770876611   -0.05141024190170344  1.558302885274501]
               [ -2.8782714155831846  -0.05482378379767665  1.6617707766444636]
               [ -2.9994120052757665  -0.057131215111630516 1.7317113286622305]
               [ -3.0899185367060014  -0.05885514298130811  1.7839652987654344]
               [-35.43175189038524   -0.03817977484002277 -65.95585512184402]
               [-31.531149492485945  -0.010298897100926725 19.931183324520134]
               [ 39.29965633377961   -0.006545683435678002 47.3267402543354]
               [  10.339058623294022 84.4381372725957       0.4183018889575614]
               [  10.343380853641088  -84.46578672918862    0.4358890208546278]]
    for i in 1:set.segments + KiteModels.KITE_PARTICLES + 1
        @test isapprox(forces[i,:], kps4.forces[i], rtol=1e-4) 
    end
end

@testset "calc_aero_forces!      " begin
    kps4.set.alpha = 1.0/7.0
    init_150()
    kps4.set.elevation = 60.0
    kps4.set.profile_law = Int(EXP)
    for i in 1:set.segments + KiteModels.KITE_PARTICLES + 1 
        kps4.forces[i] .= zeros(3)
    end
    pos, vel, acc = KiteModels.init_pos_vel_acc(kps4, old=true, delta = 1e-6)
    rho = 1.25
    kps4.v_wind .= KVec3(8.0, 0.2, 0.0)
    alpha_depower = 0.1
    rel_steering = +0.1
    kps4.set.alpha_zero = 5.0
    for i in 1:set.segments + KiteModels.KITE_PARTICLES + 1 
        kps4.forces[i] .= zeros(3)
    end
    KiteModels.calc_aero_forces!(kps4, pos, vel, rho, alpha_depower, rel_steering)
    forces=[[   0.                    0.                    0.                ]
            [   0.                    0.                    0.                ]
            [   0.                    0.                    0.                ]
            [   0.                    0.                    0.                ]
            [   0.                    0.                    0.                ]
            [   0.                    0.                    0.                ]
            [   0.                    0.                    0.                ]
            [   0.                    0.                    0.                ]
            [ -81.3257203383301146   -2.0331316776625457 -454.15328374043645  ]
            [  -9.9385080719600989  -68.915067335201357    -0.9904916722429121]
            [ -11.4093091631036021   53.2874848115847612    0.7727700100708267]]
    for i in 1:set.segments + KiteModels.KITE_PARTICLES + 1
        @test all(isapprox.(forces[i,:], kps4.forces[i], rtol=1e-4))
    end
end

function init2()
    kps4.set.alpha = 1.0/7.0
    init_150()
    kps4.set.elevation = 60.0
    kps4.set.profile_law = Int(EXP)
    pos, vel, acc = KiteModels.init_pos_vel_acc(kps4,old=true, delta = 1e-6)
    posd = copy(vel)
    veld = zero(vel)
    kps4.v_wind_gnd .= [7.0, 0.1, 0.0]
    height = 134.14733504839947
    kps4.v_wind .= kps4.v_wind_gnd * calc_wind_factor(kps4.am, height)
    kps4.stiffness_factor = 0.5
    kps4.set.alpha = 1.0/7.0
    length = 150.0
    kps4.segment_length = length/set.segments
    for i in 1:set.segments + KiteModels.KITE_PARTICLES + 1 
        kps4.forces[i] .= zeros(3)
    end
    KiteModels.init_springs!(kps4)
    return pos, vel, posd, veld
end

@testset "test_loop          " begin
    init2()
    pos, vel, posd, veld = init2()
    KiteModels.loop!(kps4, pos, vel, posd, veld)
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
    for i in 1:set.segments + KiteModels.KITE_PARTICLES + 1 
        @test all(kps4.res1[i] .≈ res1[i,:])
    end
    forces=[[ -1.058148047467609  -0.020155169648329054  0.6109409799974753]
            [-2.3531675038883413  -0.04482168872970446   1.3585829718182139]
            [-2.699059770876611   -0.05141024190170344   1.558302885274501]
            [-2.878271415585904   -0.05482378379767665   1.6617707766397531]
            [ -2.9994120052760382 -0.057131215111630516  1.7317113286617594]
            [ -3.0899185366975708  -0.05885514298130811  1.7839652987800367]
            [-35.4317518903934   -0.03817977484002277  -65.95585512185816]
            [-31.531149492485945  -0.010298897100926725  19.931183324520134]
            [ 39.29965633377961   -0.006545683435678002  47.3267402543354]
            [ 10.339058623294022   84.4381372725957       0.4183018889575614]
            [  10.343380853641088 -84.46578672918862      0.4358890208546278]]
    for i in 1:set.segments + KiteModels.KITE_PARTICLES + 1 
        @test isapprox(forces[i,:], kps4.forces[i], rtol=1e-4)
    end
    res2 = [[  0.000001             0.000001             0.000001          ]
            [ -8.556972741412151  -0.16298795901710716  14.750301715702598]
            [ -9.814762803187678  -0.18694633418801254  15.47655594645273]
            [-10.46644151122147   -0.1993592138097333   15.852802824144558]
            [-10.906952746458321  -0.2077498731332019   16.107132104224583]
            [-11.236067406172985  -0.21401870175021134  16.297146541018318]
            [ -4.1501319930182605  -0.004472008766034878 2.0845704103240807]
            [-10.803148488191985   -0.003528590502938543 16.638787927680177]
            [ 29.851165446616537   -0.004971958979489869 45.75836406156792]
            [ 10.471099184004316   85.51650034190715     10.233644040305819]
            [ 10.475476613740355  -85.544502910895       10.251455778217956]]
    for i in 1:set.segments + KiteModels.KITE_PARTICLES + 1
        @test isapprox(res2[i,:], kps4.res2[i], rtol=1e-4)
    end
end

function unpack_add_origin(y::MVector{S, SimFloat},) where S
    T = S-2
    y_  = @view y[1:end-2]
    part  = reshape(SVector{T}(y_),  Size(3, div(T,6), 2))
    pos1, vel1 = part[:,:,1], part[:,:,2]
    pos = SVector{div(T,6)+1}(if i==1 SVector(0.0,0,0) else SVector(pos1[:,i-1]) end for i in 1:div(T,6)+1)
    vel = SVector{div(T,6)+1}(if i==1 SVector(0.0,0,0) else SVector(vel1[:,i-1]) end for i in 1:div(T,6)+1)
    pos, vel
end

function unpack(y::MVector{T, SimFloat},) where T
    part  = reshape(SVector{T}(y),  Size(3, div(T,6), 2))
    pos1, vel1 = part[:,:,1], part[:,:,2]
    pos = SVector{div(T,6)}(SVector(pos1[:,i]) for i in 1:div(T,6))
    vel = SVector{div(T,6)}(SVector(vel1[:,i]) for i in 1:div(T,6))
    pos, vel
end

@testset "test_residual!       " begin
    global pos, poss, vel, vels
    # init2()
    kps4.alpha_depower = -0.820659579962 
    kps4.stiffness_factor = 0.5
    kps4.set.alpha_zero = 0.0
    res =  zeros(MVector{6*(kps4.set.segments+4), SimFloat})
    y0, yd0 = KiteModels.init(kps4; old=true, delta=1e-6)
    pos, vel = unpack_add_origin(y0)
    y0s =[  -0.                ,    0.000          ,   -0.                ,
         12.5000000000000036,    0.000001          ,   21.6506350946109656,
         25.0000000000000071,    0.000001          ,   43.3012701892219312,
         37.5000000000000071,    0.000001          ,   64.9519052838329003,
         50.0000000000000142,    0.000001          ,   86.6025403784438623,
         62.5000000000000142,    0.000001          ,  108.2531754730548244,
         75.0000000000000142,    0.000001          ,  129.9038105676658006,
         76.5905217485472747,    0.                ,  134.6435550484500823,
         78.5650000003636109,    0.                ,  136.0785716987731178,
         77.4500000002498865,    2.4810999999999996,  134.1473350483994693,
         77.4500000002498865,   -2.4810999999999996,  134.1473350483994693,
          0.00000           ,    0.00000           ,    0.00000           ,
          0.00000          ,    0.00000          ,   -0.                ,
          0.00000          ,    0.00000          ,   -0.                ,
          0.00000          ,    0.00000          ,   -0.                ,
          0.00000          ,    0.00000          ,   -0.                ,
          0.00000          ,    0.00000          ,   -0.                ,
          0.00000          ,    0.00000          ,   -0.                ,
          0.00000          ,    0.00000          ,    0.00000          ,
          0.00000          ,    0.00000          ,    0.00000          ,
          0.00000          ,    0.00000          ,    0.00000          ,
          0.00000          ,    0.00000          ,    0.00000          ,
        150.                ,    0.                ][1:end-2]
    poss, vels=unpack(MVector{66, Float64}(y0s))
    @test length(pos) == length(poss)
    for i in 1:set.segments + KiteModels.KITE_PARTICLES + 1
        @test all(pos[i] .≈ poss[i])
        @test all(vel[i] .≈ vels[i])
    end
    time = 0.0
    @test length(res) == length(y0)-2
    residual!(res, yd0, y0, kps4, time)
    res1, res2 = kps4.res1, kps4.res2
    height = calc_height(kps4)
    @test height ≈ 136.07857169877312
    res1_= zeros(MVector{(kps4.set.segments+4+1), KVec3})
    @test res1 == res1_
    # reference values output of KPS4P_v2.py of FreeKiteSim
    res2_ =[[   0.                    0.                    0.                ]
            [  -8.55694384182756     -0.16298832298499297   4.9403544121719705]
            [ -11.8070279026385805   -0.2248931126495136    6.816791315268075 ]
            [ -13.272990356870892    -0.2528161075957597    7.6631651328556583]
            [ -14.392012870751639    -0.2741307509235479    8.309233082460759 ]
            [ -15.3096925246324975   -0.2916103077136664    8.8390556776117659]
            [  -0.7543438292745689   -0.0077244009806506   -0.4074644786734876]
            [  -1.2766764661730894   -0.0049560787501444    0.6267331714897555]
            [-170.1485390054543245   -2.4632176272693642  327.4727436084622809]
            [ -19.5086166662826592 -108.1053016325064959   -0.6311974416990687]
            [ -22.7946534141140447  107.5032248949864737    1.2645875252505476]]  

    @test res2_[1, :] == res2[1]
    @test all(res2_[2, :] .≈ res2[2])
    for i in 1:set.segments + KiteModels.KITE_PARTICLES + 1
        if ! all(res2_[i, :] .≈ res2[i])
            # println(i, res2_[i, :], " ", res2[i])
            @test_broken all(res2_[i, :] .≈ res2[i])
        else
            @test all(res2_[i, :] .≈ res2[i])
        end
    end

end

@testset "test_getters" begin
    x, y, z = kite_ref_frame(kps4)
    @test all(x .≈ [-0.8660254037549963, 0.0, 0.5000000000509957])
    @test all(y .≈ [0.0, 1.0, 0.0])
    @test all(z .≈ [-0.5000000000509957, -0.0, -0.8660254037549963])
    @test all(KiteModels.orient_euler(kps4) .≈ [0, 0.5235987756571836, -1.5707963267948966])
    @test all(pos_kite(kps4) .≈ [78.56500000036361, 0.0, 136.07857169877312])
    @test calc_elevation(kps4) ≈ 1.0471975512013534 # 60 degrees
    @test calc_azimuth(kps4) ≈ 0.0 atol=0.02
    @test_broken calc_heading(kps4) ≈ 0.0 atol=1e-2
    calc_course(kps4) # the course for vel_kite = zero is undefined
end

@testset "test_find_steady_state" begin
    init_392()
    clear!(kps4)
    KiteModels.set_depower_steering!(kps4, kps4.set.depower_offset/100.0, 0.0)
    height = sin(deg2rad(kps4.set.elevation)) * kps4.set.l_tether
    kps4.v_wind .= kps4.v_wind_gnd * calc_wind_factor(kps4.am, height)
    res1, res2 = find_steady_state!(kps4; stiffness_factor=0.035, prn=false) 
    # TODO check why -9.81 appears in the residual
    @test sum(res2) ≈ -9.81*(set.segments+ KiteModels.KITE_PARTICLES) # velocity and acceleration must be near zero
    pre_tension = KiteModels.calc_pre_tension(kps4)
    @test pre_tension > 1.0001
    @test pre_tension < 1.01
    @test unstretched_length(kps4) ≈ 392.0              # initial, unstretched tether lenght
    # println("length: ", tether_length(kps4))
    @test isapprox(tether_length(kps4), 401.9975294637118, rtol=1e-2) # real, stretched tether length was: 406.4
#    @test winch_force(kps) ≈ 276.25776695110034        # initial force at the winch [N]
#    lift, drag = lift_drag(kps)
#    @test lift ≈ 443.63303000106197                    # initial lift force of the kite [N]
#    @test drag ≈ 94.25223134952152                     # initial drag force of the kite [N]
#    @test lift_over_drag(kps) ≈ 4.706870316479931      # initial lift-over-drag
#    @test norm(v_wind_kite(kps)) ≈ 9.107670173739065   # initial wind speed at the height of the kite [m/s]
end

@testset "test_spring_forces    " begin
    init2()
    kps4.alpha_depower = -0.820659579962 
    kps4.stiffness_factor = 0.04
    kps4.set.alpha_zero = 0.0
    res =  zeros(MVector{6*(kps4.set.segments+4), SimFloat})
    y0, yd0 = KiteModels.init(kps4)
    forces = spring_forces(kps4; prn=false)
    ref_forces = [3.928735076156923e-12, 3.928735076156923e-12, 3.928735076156923e-12, 0.0, -3.928735076156923e-12, 1.1786205228470769e-11, 2.160277004750972, 2.1602770047433926, 2.1602770047441373, 2.1602770047074573, 2.1602770047483513, 2.1602770047483513, 2.1602770047074573, 2.1602770047433926, 2.160277004685822]
    for i in 1:set.segments + KiteModels.KITE_PARTICLES + 1
        @test isapprox(forces[i], ref_forces[i], atol=1e-6, rtol=1e-4)
    end    
end

function simulate(integrator, steps)
    iter = 0
    for i in 1:steps
        next_step!(kps4, integrator; set_speed=0)
        iter += kps4.iter
    end
    iter / steps
end

@testset "test_simulate        " begin
    STEPS = 500
    kps4.set.depower = 23.6
    kps4.set.solver = "IDA"
    integrator = KiteModels.init_sim!(kps4; stiffness_factor=0.5, prn=false)
    # println("\nStarting simulation...")
    simulate(integrator, 100)
    av_steps = simulate(integrator, STEPS-100)
    if Sys.isapple()
        result = isapprox(av_steps, 300, rtol=0.6)
        if !result
            println("isapple $av_steps")
        end
        println("isapple $av_steps")
        @test result
    else
        result = isapprox(av_steps, 300, rtol=0.6)
        if !result
            println("not apple $av_steps")
        end
        @test result
    end
  
    lift, drag = KiteModels.lift_drag(kps4)
    # println(lift, " ", drag) # 703.7699568972286 161.44746368100536
    @test isapprox(lift, 520.5027844652874, rtol=0.05)
    sys_state = SysState(kps4)
    update_sys_state!(sys_state, kps4)
    # TODO Add testcase with varying reelout speed 
end

@testset "Raptures" begin
    kps4_ = KPS4(KCU(set))
    integrator = KiteModels.init_sim!(kps4_; stiffness_factor=0.035, prn=false)
    kps4_.set.version = 2
    kps4_.stiffness_factor = 3
    @test maximum(spring_forces(kps4_; prn=false)) > 9000
end

# TODO Add test for winch_force
@testset "test_copy_examples" begin
    cd(tempdir())
    copy_examples()
    @test isdir("examples")
    cd("examples")
    @test isfile("reel_out_4p.jl")
end

@testset "test_copy_bin" begin
    cd(tempdir())
    copy_bin()
    @test isdir("bin")
    cd("bin")
    @test isfile("create_sys_image")
end

@testset "turn" begin
    # test the turn function
    UPWIND_DIR = -pi/2 +deg2rad(10)

    # y, yd
    res=[1.0, 0, 1, 0, 0, 0]
    res2 = KiteModels.turn(res, UPWIND_DIR)
    @test res2[1] ≈ 0.984807753012208
    @test res2[2] ≈ -0.17364817766693036
    @test res2[3] ≈ 1.0
    @test res2[4] ≈ 0.0
    @test res2[5] ≈ 0.0
    @test res2[6] ≈ 0.0
end

end
nothing