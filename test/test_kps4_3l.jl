using Test, BenchmarkTools, StaticArrays, LinearAlgebra, KiteUtils
using KiteModels, KitePodModels

if ! @isdefined kcu
    const kcu = KCU(se("system_3l.yaml"))
end
if ! @isdefined kps4_3l
    kcu.set.winch_model = "AsyncMachine"
    const kps4_3l = KPS4_3L(kcu)
end

pos, vel = nothing, nothing

@testset verbose = true "KPS4_3L tests...." begin

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
    kps4_3l.set.aero_surfaces = 3
    
    kps4_3l.set.sim_settings = "3l_settings.yaml"
    kps4_3l.set.sim_time = 100.0
    kps4_3l.set.abs_tol = 0.006
    kps4_3l.set.rel_tol = 0.01
    kps4_3l.set.max_iter = 10000
    kps4_3l.set.physical_model = "KPS4_3L"
    kps4_3l.set.version = 2
    kps4_3l.set.cl_list = [0.0, 0.5, 0.0, 0.08, 0.125, 0.15, 0.0, 1.0, 1.0, 0.0, -0.5, 0.0]
    kps4_3l.set.cd_list = [0.5, 0.5, 0.5,  1.0,   0.2,  0.1, 0.2, 1.0, 0.5, 0.5,  0.5]
    kps4_3l.set.alpha_zero = 10.0
    kps4_3l.set.d_tether = 1.0
    kps4_3l.set.v_wind_ref = [15.51, 0.0]
    kps4_3l.set.depower = 25.0
    kps4_3l.set.alpha = 0.08163
    KiteModels.clear!(kps4_3l)
    # kps4_3l.set.
end

function init_100()
    KiteModels.clear!(kps4_3l)
    kps4_3l.set.l_tether = 100.0
    kps4_3l.set.elevation = 70.0
    kps4_3l.set.v_wind = 15.51
    kps4_3l.set.mass = 0.5
    kps4_3l.set.c_s = 0.6
    kps4_3l.set.damping = 473.0     # unit damping coefficient
    kps4_3l.set.c_spring = 614600.0 # unit spring coefficent
    kps4_3l.set.width = 3.0
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

function init3()
    KiteModels.clear!(kps4_3l)
    kps4_3l.set.l_tether = 150.0 # - kps4_3l.set.height_k - kps4_3l.set.h_bridle
    kps4_3l.set.area = 10.18
    kps4_3l.set.mass = 6.21
    kps4_3l.set.c_s = 0.6
    kps4_3l.set.damping = 473.0     # unit damping coefficient
    kps4_3l.set.c_spring = 614600.0 # unit spring coefficent
    kps4_3l.set.width = 4.9622
    kps4_3l.set.elevation = 70.7 
    kps4_3l.set.profile_law = Int(EXPLOG)
    pos, vel = KiteModels.init_pos_vel(kps4_3l)
    posd = copy(vel)
    veld = zero(vel)
    height = 134.14733504839947
    kps4_3l.v_wind .= kps4_3l.v_wind_gnd * calc_wind_factor(kps4_3l.am, height)
    kps4_3l.stiffness_factor = 1.0
    KiteModels.init_springs!(kps4_3l)
    return pos, vel, posd, veld
end

set_defaults()

@testset "calc_rho              " begin
    @test isapprox(calc_rho(kps4_3l.am, 0.0), 1.225, atol=1e-5) 
    @test isapprox(calc_rho(kps4_3l.am, 100.0), 1.210756, atol=1e-5) 
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
    init_50()
    particles = KiteModels.get_particles(kps4_3l.set.width, kps4_3l.set.radius, kps4_3l.set.middle_length, kps4_3l.set.tip_length, kps4_3l.set.bridle_center_distance)
    @test particles[1] == [75.0, 0.0, 129.90381057] # E
    @test particles[2] == [76.95106065350437, 0.6180085233533619, 133.28314675005853]   # C
    @test particles[3] == [76.95106065350437, -0.6180085233533616, 133.28314675005853]   # D
    @test particles[4] == [77.64618475974454, 1.1102230246251565e-16, 132.881816660146]   # A
end

@testset "init_springs!          " begin
    init_50()
    sp = KiteModels.init_springs!(kps4_3l)
    # test springs
    @test length(sp) == 6*3 + 6
    @test sp[1].p1 == 1
    @test sp[1].p2 == 4
    for i in 1:6*3
        @test sp[i].length   ≈ 8.333333333333334
        @test sp[i].c_spring ≈ 5183.627878423158
        @test sp[i].damping  ≈ 56.76
    end
    @test sp[18].p1 == 18
    @test sp[18].p2 == 21
    @test sp[19].p1 == 21
    @test sp[19].p2 == 24
    @test sp[19].length ≈ 3.9830222651726808
    @test sp[19].c_spring ≈ 67782.8544993504
    @test sp[19].damping ≈ 118.75404366575728
    @test sp[20].length ≈ 3.949967399446711
    @test sp[21].length ≈ 3.949967399446711
    @test sp[22].length ≈ 1.2357698432973823
    @test sp[23].length ≈ 1.0128116611547708
    @test sp[24].c_spring ≈ 266565.4721629595
end

@testset "init_masses!           " begin
    init_50()
    m = KiteModels.init_masses!(kps4_3l)
    @test m[1] ≈ 0.0
    for i in 4:18
        @test m[i] ≈ 0.004738568919164605
    end
    @test m[19] ≈ 0.0023692844595823025
    @test m[20] ≈ 0.0023692844595823025
    @test m[21] ≈ 0.01186537657358817
    @test m[22] ≈ 0.125
    @test m[23] ≈ 0.125
    @test m[end] ≈ 0.25
end

@testset "calc_particle_forces!  " begin
    init_50()
    pos1 = KVec3(1.0, 2.0, 3.0)
    pos2 = KVec3(2.0, 3.0, 4.0)
    vel1 = KVec3(3.0, 4.0, 5.0)
    vel2 = KVec3(4.0, 5.0, 6.0)
    rho = kps4_3l.set.rho_0
    for i in 1:kps4_3l.num_A 
        kps4_3l.forces[i] .= zeros(3)
    end
    bytes = 0
    for i in 1:length(kps4_3l.springs)
        spring = kps4_3l.springs[i]
        kps4_3l.stiffness_factor = 0.5
        kps4_3l.v_wind_tether .= KVec3(8.0, 0.1, 0.0)
        bytes = @allocated KiteModels.calc_particle_forces!(kps4_3l, pos1, pos2, vel1, vel2, spring, se().d_tether/1000.0, rho, i)
    end
    # @test bytes == 0
    res=[[-931.0208419755234   -931.0559150592786   -931.0602499347989]
        [-931.0208419755234   -931.0559150592786   -931.0602499347989]
        [-931.0208419755234   -931.0559150592786   -931.0602499347989]
        [0.04965402868720048   -0.020492138823328787   -0.029161889863871693]
        [0.04965402868720048   -0.020492138823328787   -0.029161889863871693]
        [0.04965402868720048   -0.020492138823328787   -0.029161889863871693]
        [0.04965402868720048   -0.020492138823328787   -0.029161889863871693]
        [0.04965402868720048   -0.020492138823328787   -0.029161889863871693]
        [0.04965402868720048   -0.020492138823328787   -0.029161889863871693]
        [0.04965402868720048   -0.020492138823328787   -0.029161889863871693]
        [0.04965402868720048   -0.020492138823328787   -0.029161889863871693]
        [0.04965402868720048   -0.020492138823328787   -0.029161889863871693]
        [0.04965402868720048   -0.020492138823328787   -0.029161889863871693]
        [0.04965402868720048   -0.020492138823328787   -0.029161889863871693]
        [0.04965402868720048   -0.020492138823328787   -0.029161889863871693]
        [0.04965402868720048   -0.020492138823328787   -0.029161889863871693]
        [0.04965402868720048   -0.020492138823328787   -0.029161889863871693]
        [0.04965402868720048   -0.020492138823328787   -0.029161889863871693]
        [931.0704960042106   931.0354229204553   931.031088044935]
        [931.0704960042106   931.0354229204553   931.031088044935]
        [-11867.4612536364   -11867.601545971422   -11867.618885473503]
        [96000.29073757448   96000.18551832321   96000.17251369666]
        [28808.982257032563   28808.8770377813   28808.86403315474]
        [-112010.44332079432   -112010.54854004558   -112010.56154467215]]
    for i in 1:kps4_3l.num_A
        @test all(res[i,:] .≈ kps4_3l.forces[i])
        # println(kps4_3l.forces[i])
    end
end

@testset "init                  " begin
    init_50()
    kps4_3l.set.elevation = 60.0
    pos, vel, acc = KiteModels.init_pos_vel_acc(kps4_3l, delta = 1e-6)
    pos1 = [[0.0   1.0e-6   0.0]
            [0.0   1.0e-6   0.0]
            [0.0   1.0e-6   0.0]
            [4.491843442217564   0.10300158722556031   7.780101061565895]
            [4.491843442217564   -0.10300158722556031   7.780101061565895]
            [4.166666666666668   1.0e-6   7.216878364870322]
            [8.983686884435128   0.20600317445112062   15.56020212313179]
            [8.983686884435128   -0.20600317445112062   15.56020212313179]
            [8.333333333333336   1.0e-6   14.433756729740644]
            [13.475530326652692   0.30900476167668095   23.340303184697685]
            [13.475530326652692   -0.30900476167668095   23.340303184697685]
            [12.500000000000004   1.0e-6   21.650635094610966]
            [17.967373768870257   0.41200634890224125   31.12040424626358]
            [17.967373768870257   -0.41200634890224125   31.12040424626358]
            [16.66666666666667   1.0e-6   28.867513459481287]
            [22.45921721108782   0.5150079361278016   38.90050530782948]
            [22.45921721108782   -0.5150079361278016   38.90050530782948]
            [20.83333333333334   1.0e-6   36.08439182435161]
            [26.951060653305383   0.6180095233533619   46.68060636939537]
            [26.951060653305383   -0.6180095233533619   46.68060636939537]
            [25.000000000000007   1.0e-6   43.30127018922193]
            [26.951060653305383   0.6180095233533619   46.68060636939537]
            [26.951060653305383   -0.6180095233533619   46.68060636939537]
            [27.64618475956919   1.000000000139778e-6   46.27927627952376]
           ]
    for i in 1:kps4_3l.num_A
        @test all(pos1[i,:] .≈ pos[i])
        # println(pos[i])
    end
end


@testset "initial_residual      " begin
    function split_res(res, num_particles)
        pos=res[1:3*particles]
        vel=res[3*particles+1:6*particles]
        pos, vel
    end
    init3()
    res = zeros(MVector{6*(kps4_3l.num_A-5)+4+6, SimFloat})
    X00 = zeros(SimFloat, 5*kps4_3l.set.segments+3)
    y0, yd0 = KiteModels.init(kps4_3l, X00)
    num_particles = div(length(y0)-6-4, 6)
    residual!(res, yd0, y0, kps4_3l, 0.0)
    res_pos, res_vel = res[1:3*num_particles+2], res[3*num_particles+3:6*num_particles+4]
    @test res_pos == zeros(length(res_pos))
    # in the second test we check if the norm of the accelerations of the tether particles is low
    # @test (norm(res_vel)) < 15.0
    # println(res_vel)
end

@testset "inner_loop!            " begin
    init_50()
    kps4_3l.set.elevation = 60.0
    kps4_3l.set.profile_law = Int(EXP)
    for i in 1:kps4_3l.num_A 
        kps4_3l.forces[i] .= zeros(3)
    end
    pos, vel, acc = KiteModels.init_pos_vel_acc(kps4_3l, delta = 1e-6)
    v_wind_gnd = KVec3(7.0, 0.1, 0.0)
    kps4_3l.stiffness_factor = 0.5
    segments = kps4_3l.set.segments
    d_tether = kps4_3l.set.d_tether/1000.0
    KiteModels.inner_loop!(kps4_3l, pos, vel, v_wind_gnd, d_tether)
    forces =  [[-691.7059635853744   -15.86212132213054   -1198.250303656443]
                [-691.7059267098825   15.865405658555337   -1198.2502791771315]
                [-719.8767233654121   0.0013636636282290742   -1247.0283971963331]
                [0.17151129393744213   0.0018017696171863662   -0.09904480166824214]
                [0.1715834423043816   0.0044245018272306424   -0.09900631367599999]
                [0.15717840676779815   0.0029938448442505183   -0.09074699545817566]
                [0.19480217915588582   0.0022213468985956553   -0.11249849922182875]
                [0.1948872556117749   0.005200264806310528   -0.11244936262801275]
                [0.17853582897430442   0.003400650929770767   -0.10307770891836299]
                [0.20846142849427451   0.0023771057706643006   -0.12038673261440636]
                [0.20855247041777147   0.005564901162090408   -0.12033415061682717]
                [0.1910670017234679   0.0036393388566165795   -0.11031258487832929]
                [0.21847280093288646   0.0024912672782040346   -0.1261683125815125]
                [0.21856821517053504   0.005832156755717577   -0.12611320531914316]
                [0.20025615129509333   0.003814369514879717   -0.11561794285694305]
                [0.2264609993015938   0.0025823581803336992   -0.13078150702381208]
                [0.22655990226348877   0.006045403330494636   -0.13072438482117832]
                [0.2075919258348904   0.003954097961792049   -0.11985325426326199]
                [691.899119571082   15.864477896223933   1198.138752653517]
                [691.8991698100529   -15.860095319097438   1198.1387816694394]
                [8610.386516090239   0.0015829749097520107   11921.484451159264]
                [-4005.3977950461667   14519.255397533867   -1104.130875838349]
                [-4005.375232614203   -14519.232357143173   -1104.1559349816287]
                [120.79342186625172   -0.0162820172326974   -8466.429543362992]]
    for i in 1:kps4_3l.num_A
        @test isapprox(forces[i,:], kps4_3l.forces[i], rtol=1e-4)
        # println(kps4_3l.forces[i])
    end
end

@testset "calc_aero_forces!      " begin
    init_50()
    kps4_3l.set.elevation = 60.0
    kps4_3l.set.profile_law = Int(EXP)
    for i in 1:kps4_3l.num_A
        kps4_3l.forces[i] .= zeros(3)
    end
    pos, vel, acc = KiteModels.init_pos_vel_acc(kps4_3l, delta = 1e-6)
    rho = 1.25
    kps4_3l.v_wind .= KVec3(8.0, 0.2, 0.0)
    for i in 1:kps4_3l.num_A
        kps4_3l.forces[i] .= zeros(3)
    end
    KiteModels.calc_aero_forces!(kps4_3l, pos, vel)
    forces=[[0.0  0.0  0.0]
            [0.0  0.0  0.0]
            [0.0  0.0  0.0]
            [0.0  0.0  0.0]
            [0.0  0.0  0.0]
            [0.0  0.0  0.0]
            [0.0  0.0  0.0]
            [0.0  0.0  0.0]
            [0.0  0.0  0.0]
            [0.0  0.0  0.0]
            [0.0  0.0  0.0]
            [0.0  0.0  0.0]
            [0.0  0.0  0.0]
            [0.0  0.0  0.0]
            [0.0  0.0  0.0]
            [0.0  0.0  0.0]
            [0.0  0.0  0.0]
            [0.0  0.0  0.0]
            [4.937589983558479  0.0  14.099543704451493]
            [4.733790218796835  0.0  13.517582930109748]
            [0.0  0.0  0.0]
            [10.940304617086838  28.40311527705575  64.64855386496264]
            [11.581670317238302  -26.584094467198202  61.340408835496966]
            [0.0  0.0  0.0]]
    for i in 1:kps4_3l.num_A
        @test all(forces[i,:] .≈ kps4_3l.forces[i])
        # println(kps4_3l.forces[i])
    end
end

function init2()
    init_50()
    kps4_3l.set.elevation = 60.0
    kps4_3l.set.profile_law = Int(EXP)
    pos, vel, acc = KiteModels.init_pos_vel_acc(kps4_3l, delta = 1e-6)
    posd = copy(vel)
    veld = zero(vel)
    kps4_3l.v_wind_gnd .= [7.0, 0.1, 0.0]
    height = 134.14733504839947
    kps4_3l.v_wind .= kps4_3l.v_wind_gnd * calc_wind_factor(kps4_3l.am, height)
    kps4_3l.stiffness_factor = 0.5
    lengths = [155.0, 155.0, 150.0]
    kps4_3l.segment_lengths .= lengths/se().segments
    for i in 1:kps4_3l.num_A 
        kps4_3l.forces[i] .= zeros(3)
    end
    KiteModels.init_springs!(kps4_3l)
    return pos, vel, posd, veld
end

@testset "test_loop          " begin
    init2()
    pos, vel, posd, veld = init2()
    KiteModels.loop!(kps4_3l, pos, vel, posd, veld)
    for i in 1:kps4_3l.num_A 
        @test all(kps4_3l.res1[i] .≈ zeros(3))
    end
    forces=[[-28813.450404866944 -660.7101060126089 -49906.54047342375]
            [-28813.45025799939 660.7259090755075 -49906.5402584326]
            [-30729.92840691776 0.0013636636282290742 -53225.9626508769]
            [0.17156628970769816 -0.004457593638107937 -0.09894954620540375]
            [0.17152844600423123 -0.0018348613912166911 -0.09910157005651854]
            [0.1571784067673434 0.0029938448442505183 -0.09074699545453768]
            [0.19480217916861875 0.0022213468986365115 -0.11249849921296118]
            [0.1948872556204151 0.005200264806262567 -0.11244936261937255]
            [0.1785358289744181 0.003400650929770767 -0.10307770891813561]
            [0.20846142847949523 0.0023771057706198917 -0.12038673261122312]
            [0.20855247039798996 0.00556490116207442 -0.12033415061887354]
            [0.19106700172415003 0.0036393388566165795 -0.11031258488219464]
            [0.21847280094880261 0.0024912672780601497 -0.12616831258492311]
            [0.2185682151830406 0.005832156755786855 -0.12611320532596437]
            [0.2002561513108958 0.003814369514879717 -0.11561794285807991]
            [0.2264609992998885 0.0025823581802342233 -0.13078150703950087]
            [0.22655990226121503 0.006045403330631416 -0.1307243848350481]
            [0.2075919257949863 0.003954097961792049 -0.11985325427667703]
            [18712.336826052975 0.0 53434.046117617705]
            [18712.33684059229 0.0 53434.046159135476]
            [30774.854649875913 0.0030858773048274557 53292.16389482203]
            [10106.584909812285 612.8007690783015 -3571.248089177784]
            [10106.691184520216 -612.6931339600861 -3571.32147151677]
            [-55.05875505245474 -0.09799707002383684 1.2534231229097088]]
    for i in 1:kps4_3l.num_A 
        @test isapprox(forces[i,:], kps4_3l.forces[i], rtol=1e-4)
        # println(kps4_3l.forces[i]')
    end
    res2 = [[0.0 0.0 0.0]
            [0.0 0.0 0.0]
            [0.0 0.0 0.0]
            [-11.679468076955835 0.3034531007520357 16.546043940241628]
            [-11.676891846349527 0.12490918280546429 16.556393046228916]
            [-11.056671993074616 -0.21060119593929524 16.193572551304943]
            [-13.261263834507382 -0.1512193930093825 17.468396253143812]
            [-13.267055460033287 -0.3540108426890771 17.46505125292724]
            [-12.559054011176961 -0.23921785865328948 17.06097321410335]
            [-14.191124627984719 -0.16182276256485825 18.005392013298923]
            [-14.197322355886486 -0.378833660065718 18.001812465695977]
            [-13.440555927578416 -0.25600829552129734 17.569908583091962]
            [-14.87265423010529 -0.16959436900361324 18.3989762008001]
            [-14.879149605768937 -0.3970272293296467 18.395224744305274]
            [-14.086964139530433 -0.26832077925840525 17.9431125374496]
            [-15.416455158556053 -0.17579543149587662 18.713021911392936]
            [-15.423188031205438 -0.41154410538764463 18.70913328639249]
            [-14.602997187287517 -0.2781499386027814 18.241044359683094]
            [-7.897885266740787e6 0.0 -2.255282007253699e7]
            [-7.897885272877371e6 0.0 -2.2552820090060327e7]
            [-2.5936686003189688e6 -0.26007411443615613 -4.491391163606531e6]
            [-80852.67927849828 -4902.406152626412 28579.794713422274]
            [-80853.52947616173 4901.545071680689 28580.381772134162]
            [220.23502020981897 0.39198828009534736 4.796307508361165]]
    for i in 1:kps4_3l.num_A
        @test isapprox(res2[i,:], kps4_3l.res2[i], rtol=1e-4)
        # println(kps4_3l.res2[i]')
    end
end

@testset "test_residual!       " begin
    
    kps4_3l.stiffness_factor = 0.04
    res =  zeros(MVector{6*(kps4_3l.num_A-5)+4+6, SimFloat})
    y0, yd0 = KiteModels.init(kps4_3l; delta=1e-6)
    y0s = [4.490372205878421, 0.10298757652176786, 7.780950480262341, 4.490372205878421, -0.10298757652176786, 7.780950480262341, 4.166666666666668, 1.0e-6, 7.216878364870322, 8.980744411756842, 0.20597515304353572, 15.561900960524682, 8.980744411756842, -0.20597515304353572, 15.561900960524682, 8.333333333333336, 1.0e-6, 14.433756729740644, 13.471116617635262, 0.3089627295653036, 23.342851440787022, 13.471116617635262, -0.3089627295653036, 23.342851440787022, 12.500000000000004, 1.0e-6, 21.650635094610966, 17.961488823513683, 0.41195030608707145, 31.123801921049363, 17.961488823513683, -0.41195030608707145, 31.123801921049363, 16.66666666666667, 1.0e-6, 28.867513459481287, 22.451861029392106, 0.5149378826088393, 38.904752401311704, 22.451861029392106, -0.5149378826088393, 38.904752401311704, 20.83333333333334, 1.0e-6, 36.08439182435161, 25.000000000000007, 1.0e-6, 43.30127018922193, 26.942233235270525, 0.6179254591306071, 46.685702881574045, 26.942233235270525, -0.6179254591306071, 46.685702881574045, 27.646090205747516, 0.013239546805801463, 46.27933087019816, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 0.0, 0.0, 0.0, 0.0, 53.90566405419164, 53.90566405419164, 50.0, 0.0, 0.0, 0.0]
    yd0s = [1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    # println(y0)
    # println(yd0)
    @test all(y0 .≈ y0s)
    @test all(yd0 .≈ yd0s)

    time = 0.0
    @test length(res) == length(y0)
    residual!(res, yd0, y0, kps4_3l, time)
    res1, res2 = kps4_3l.res1, kps4_3l.res2
    height = calc_height(kps4_3l)
    @test height ≈ 46.27933087019816
    res1_= zeros(MVector{(kps4_3l.num_A), KVec3})
    @test res1 == res1_
    res2_ =[[0.0   0.0   0.0]
            [0.0  0.0  0.0]
            [0.0  0.0  0.0]
            [-33.576530377010556  -0.38282988815729174  29.208411662359655]
            [-33.591305552166176  -0.8963637655564114  29.199734330656238]
            [-33.16183456642719  -0.6318028861961738  28.97488975075343]
            [-38.14451435598529  -0.43500453818836243  31.82888830200402]
            [-38.16116195058352  -1.017913319895357  31.8192649533944]
            [-37.67716103350859  -0.7176525759598683  31.56292064234734]
            [-40.81914902429016  -0.46550669426878305  33.37281877027291]
            [-40.83696392819467  -1.0892881485455548  33.3625206441397]
            [-40.32166678257849  -0.7680238865638919  33.08972674844635]
            [-42.779486884607266  -0.4878628493400021  34.5044220444638]
            [-42.79815735246308  -1.1416013847855135  34.493629349052206]
            [-42.260891416297945  -0.8049613377752156  34.20933861018946]
            [-44.34366362984511  -0.5057010931465664  35.407341672896614]
            [-44.36301676237147  -1.1833427378553365  35.39615435476668]
            [-43.80899056925542  -0.8344488158083441  35.10313407824504]
            [-345.2546735104606  0.00017776169578435534  -601.6224946605902]
            [-331.46918359703903  0.00017066394374142083  -577.6006306044044]
            [-320.38753260380525  -0.9850270455038719  -425.82794106931027]
            [-254.62802485602307  -315.80551879718183  -801.6512583230729]
            [-1190.464532876995  -528.2876381560295  -250.28643973323926]
            [484.9562144240868  419.3822884418879  -260.4334085557314]]  

    # for r in res2
    #     println(r)
    # end
    @test res2_[1, :] == res2[1]
    @test all(res2_[4, :] .≈ res2[4])
    for i in 1:kps4_3l.num_A
        if ! all(res2_[i, :] .≈ res2[i])
            # println(i, res2_[i, :], " ", res2[i])
            @test_broken all(res2_[i, :] .≈ res2[i])
        else
            @test all(res2_[i, :] .≈ res2[i])
        end
        # println(res2[i]')
    end

end

@testset "test_getters" begin
    x, y, z = kite_ref_frame(kps4_3l)
    @test all(x .≈ [-0.8673285322802265, 0.0, 0.49773609181228007])
    @test all(y .≈ [0.0, 1.0, 0.0])
    @test all(z .≈ [-0.49773609181228007, 2.5626999001639186e-7, -0.8673285322802265])
    @test all(orient_euler(kps4_3l) .≈ [1.5707963267948966, -0.5209866063773279, 1.5707963267948966])
    @test all(pos_kite(kps4_3l) .≈ [27.646090205747516, 0.013239546805801463, 46.27933087019816])
    @test calc_elevation(kps4_3l) ≈ 1.0323095578497548
    @test calc_azimuth(kps4_3l) ≈ -0.000478893966385608
    @test calc_heading(kps4_3l) ≈ 0.0004154220077356996
    calc_course(kps4_3l) # the course for vel_kite = zero is undefined
end

@testset "test_find_steady_state" begin
    init_50()
    clear!(kps4_3l)
    height = sin(deg2rad(kps4_3l.set.elevation)) * kps4_3l.set.l_tether
    kps4_3l.v_wind .= kps4_3l.v_wind_gnd * calc_wind_factor(kps4_3l.am, height)
    res1, res2 = find_steady_state!(kps4_3l; stiffness_factor=0.035, prn=true) 
    @test sum(res2) ≈ 0.0 # not moving
    pre_tension = KiteModels.calc_pre_tension(kps4_3l)
    @test pre_tension > 1.000001
    @test pre_tension < 1.01
    @test unstretched_length(kps4_3l) ≈ 50.0              # initial, unstreched tether lenght
    @test isapprox(tether_length(kps4_3l), 51.31951395939359, rtol=1e-2)
#    @test winch_force(kps) ≈ 276.25776695110034        # initial force at the winch [N]
#    lift, drag = lift_drag(kps)
#    @test lift ≈ 443.63303000106197                    # initial lift force of the kite [N]
#    @test drag ≈ 94.25223134952152                     # initial drag force of the kite [N]
#    @test lift_over_drag(kps) ≈ 4.706870316479931      # initlal lift-over-drag
#    @test norm(v_wind_kite(kps)) ≈ 9.107670173739065   # inital wind speed at the height of the kite [m/s]
end

@testset "test_spring_forces    " begin
    init2()
    kps4_3l.stiffness_factor = 0.04
    res =  zeros(MVector{6*(kps4_3l.set.segments+4), SimFloat})
    y0, yd0 = KiteModels.init(kps4_3l)
    forces = spring_forces(kps4_3l)
    ref_forces = [-3493.569454791926, -3493.569454791926, -3455.751918948772, -3493.569454791926, -3493.569454791926, -3455.751918948772, -3493.569454791926, -3493.569454791926, -3455.751918948772, -3493.569454791926, -3493.569454791926, -3455.751918948772, -3493.5694547919256, -3493.5694547919256, -3455.751918948772, -3493.5694547919265, -3493.5694547919265, -3455.751918948773, 2.1602770047500774, 2.1602770047559905, 2.1602770047559905, 2.160277004745103, 2.160277004793414, 2.160277004793414]
    # println("forces ", forces)
    for i in 1:kps4_3l.num_A
        @test isapprox(forces[i], ref_forces[i], atol=1e-6, rtol=1e-4)
    end    
end

function simulate(integrator, steps)
    iter = 0
    for i in 1:steps
        KiteModels.next_step!(kps4_3l, integrator, set_values=[0.0, 0.0, 0.0])
        iter += kps4_3l.iter
    end
    return iter/steps
end

@testset "test_simulate        " begin
    STEPS = 50

    kps4_3l.set.solver = "DFBDF"
    # println("finding steady state")
    init_50()
    integrator = KiteModels.init_sim!(kps4_3l; stiffness_factor=0.035, prn=false)
    # println("\nStarting simulation...")
    simulate(integrator, STEPS)
    av_steps = simulate(integrator, STEPS-10)
    if Sys.isapple()
        println("isapple $av_steps")
        @test isapprox(av_steps, 835.25, rtol=0.6)
    else
        println("not apple $av_steps")
        @test isapprox(av_steps, 835.25, rtol=0.6)
    end
  
    lift, drag = KiteModels.lift_drag(kps4_3l)
    # println(lift, " ", drag) # 703.7699568972286 161.44746368100536
    @test isapprox(lift, 404.2596735903995, rtol=0.05)
    sys_state = SysState(kps4_3l)
    update_sys_state!(sys_state, kps4_3l)
    # TODO Add testcase with varying reelout speed 
end

@testset "Steady state history" begin
    STEPS = 50
    kps4_3l.set.solver = "DFBDF"
    # println("finding steady state")
    init_50()
    steady_state_history = load_history()

    integrator = KiteModels.init_sim!(kps4_3l; stiffness_factor=0.035, prn=false, steady_state_history=steady_state_history)
    # println("\nStarting simulation...")
    simulate(integrator, STEPS)
    lift1, drag1 = KiteModels.lift_drag(kps4_3l)
    save_history(steady_state_history)

    steady_state_history2 = load_history()
    init_50()
    integrator = KiteModels.init_sim!(kps4_3l; stiffness_factor=0.035, prn=false, steady_state_history=steady_state_history2)
    simulate(integrator, STEPS)
    lift2, drag2 = KiteModels.lift_drag(kps4_3l)

    @test lift2 == lift1
    @test drag2 == drag1
  
    @test isapprox(lift1, 849.1874782889682, rtol=0.05)
end

# @testset "Raptures" begin
#     kps4_3l_ = KPS4_3L(KCU(se()))
#     integrator = KiteModels.init_sim!(kps4_3l_; stiffness_factor=0.035, prn=false)
#     kps4_3l_.set.version = 2
#     kps4_3l_.stiffness_factor = 3
#     @test maximum(spring_forces(kps4_3l_)) > 20000
# end

# # TODO Add test for winch_force
# @testset "test_copy_examples" begin
#     cd(tempdir())
#     copy_examples()
#     @test isdir("examples")
#     cd("examples")
#     @test isfile("reel_out_4p.jl")
# end

# @testset "test_copy_bin" begin
#     cd(tempdir())
#     copy_bin()
#     @test isdir("bin")
#     cd("bin")
#     @test isfile("create_sys_image")
# end

# println(kps4_3l.set)

end
nothing