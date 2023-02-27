using Test, BenchmarkTools, StaticArrays, LinearAlgebra, KiteUtils

using KiteModels, KitePodModels

include("../src/consts.jl")

const SEGMENTS = se().segments
if ! @isdefined kcu
    const kcu = KCU(se())
    const kps = KPS3(kcu)
end
res1 = zeros(SVector{SEGMENTS+1, KiteModels.KVec3})
res2 = deepcopy(res1)
if ! @isdefined res3
    if USE_WINCH
        const res3 = vcat(reduce(vcat, vcat(res1, res2)), zeros(2))
    else
        const res3 = vcat(reduce(vcat, vcat(res1, res2)))
    end
end

@testset verbose = true "KPS3 tests...." begin

function set_defaults()
    KiteModels.clear!(kps)
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
    KiteModels.clear!(kps)
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
    @test isapprox(calc_rho(kps.am, 0.0), 1.225, atol=1e-5) 
    @test isapprox(calc_rho(kps.am, 100.0), 1.210756, atol=1e-5) 
end

@testset "calc_wind_factor     " begin
    set_defaults()
    @test isapprox(calc_wind_factor(kps.am, 6.0, Int(EXP)),   1.0, atol=1e-5) 
    @test isapprox(calc_wind_factor(kps.am, 10.0, Int(EXP)),  1.0757037, atol=1e-5) 
    @test isapprox(calc_wind_factor(kps.am, 100.0, Int(EXP)), 1.494685, atol=1e-5)
end
 
@testset "calc_cl              " begin
    @test isapprox(kps.calc_cl(-5.0), 0.150002588978, atol=1e-4) 
    @test isapprox(kps.calc_cl( 0.0), 0.200085035326, atol=1e-4) 
    @test isapprox(kps.calc_cl(10.0), 0.574103590856, atol=1e-4)
    @test isapprox(kps.calc_cl(20.0), 1.0, atol=1e-4)
end

@testset "test_calc_drag       " begin
    v_segment = KVec3(1.0, 2, 3)
    unit_vector = KVec3(2.0, 3.0, 4.0)
    rho = SimFloat(calc_rho(kps.am, 10.0))
    last_tether_drag = KVec3(0.0, 0.0, 0.0)
    v_app_perp = KVec3(0, -3.0, -4.0)
    area = 20.0   
    kps.v_wind_tether .= [0.1, 0.2, 0.3]
    v_app_norm = KiteModels.calc_drag(kps, v_segment, unit_vector, rho, v_app_perp, area)
    @test v_app_norm ≈ 3.3674916481
    @test kps.last_tether_drag ≈ [-38506.7140169,  -57266.39520463, -76026.07639235]
    @test v_app_perp ≈ [ 35.1, 52.2, 69.3]
end

@testset "test_calc_aero_forces" begin
    set_defaults()
    kps.v_apparent .= KVec3(35.1, 52.2, 69.3)
    kps.v_wind .= kps.v_wind_gnd
    pos_kite = KVec3(30.0, 5.0, 100.0)
    v_kite = KVec3(3.0, 5.0, 2.0)
    rho = SimFloat(calc_rho(kps.am, 10.0))
    rel_steering = 0.1
    kps.beta = 0.1
    kps.psi = 0.2
    kps.param_cl = 0.2
    kps.param_cd = 1.0
    KiteModels.calc_aero_forces(kps, pos_kite, v_kite, rho, rel_steering)
    @test kps.v_apparent ≈ [5.0,  -5, -2]                                # same with Python
    @test kps.kite_y ≈ [ 0.64101597,  0.73258967, -0.22893427]           # same with Python
    @test kps.cor_steering ≈ 0.0250173783309
    @test kps.steering_force ≈ [-15.88482337, -18.15408385, 5.6731512 ]
    @test kps.last_force ≈ [-555.24319976, 544.82004621, 80.49946362]    # same with Python
end

@testset "test_calc_res        " begin
    set_defaults()
    i = 2
    pos1 = KVec3(30.0, 5.0, 100.0)
    pos2 = KVec3(30.0+10, 5.0+11, 100.0+20)
    vel1 = KVec3(3.0, 5.0, 2.0)
    vel2 = KVec3(3.0+0.1, 5.0+0.2, 2.0+0.3)
    mass = 9.0
    veld = KVec3(0.1, 0.3, 0.4)
    result = KVec3(0, 0, 0)
    kps.c_spring = 0.011
    kps.damping = 0.01
    kps.last_tether_drag = KVec3(5.0,6,7)
    kps.last_force = KVec3(-1.0, -2, -3)
    kps.v_app_perp = KVec3(0.1,0.22,0.33)
    kps.v_wind_tether .= [0.1, 0.2, 0.3]
    kps.segment_length = 10.0
    KiteModels.calc_res(kps, pos1, pos2, vel1, vel2, mass, veld, result, i)
    @test result ≈ [ -6.99178740e-03, 1.01297086e-01, 9.83843908e+00]
    i = SEGMENTS+1
    KiteModels.calc_res(kps, pos1, pos2, vel1, vel2, mass, veld, result, i)
    @test result ≈ [0.14263997, 0.41669208, 10.12449937]
end

@testset "test_calc_loop       " begin
    set_defaults()
    kps.last_tether_drag = KVec3(5.0,6,7)
    kps.last_force = KVec3(-1.0, -2, -3)
    kps.v_app_perp = KVec3(0.1,0.22,0.33)
    kps.v_wind_tether .= [0.1, 0.2, 0.3]
    kps.segment_length = 10.0
    kps.c_spring = kps.set.c_spring / kps.segment_length
    kps.damping  = kps.set.damping / kps.segment_length
    pos  = zeros(SVector{SEGMENTS+1, KVec3})
    for i in 1:SEGMENTS+1
        pos[i][3] = 5.0 * (i-1)
    end
    vel  = zeros(SVector{SEGMENTS+1, KVec3})
    posd = zeros(SVector{SEGMENTS+1, KVec3})
    veld = zeros(SVector{SEGMENTS+1, KVec3})
    res1 = zeros(SVector{SEGMENTS+1, KVec3})
    res2 = zeros(SVector{SEGMENTS+1, KVec3})
    @test kps.c_spring ≈ 61460.0
    @test kps.damping  ≈    94.6
    KiteModels.loop(kps, pos, vel, posd, veld, res1, res2)
    @test sum(res1) ≈ [0.0, 0.0, 0.0]
    @test isapprox(res2[7], [ -5.02874357e-02, -1.00574871e-01, -7.62063430e+02])
    @test isapprox(res2[6], [ -6.15356097e-03, -1.23071219e-02, 9.81000000e+00]) 
    @test isapprox(res2[5], [ -2.38000593e-03, -4.76001187e-03, 9.81000000e+00]) 
    @test isapprox(res2[2], [ -2.38139816e-03, -4.76279631e-03, 9.81000000e+00], rtol=1e-5) 
    @test isapprox(res2[1], [0.0,0.0,0.0], rtol=1e-5)
end

@testset "test_calc_alpha      " begin
    v_app = KVec3(10,2,3)
    vec_z = normalize(KVec3(3,2,0))
    alpha = KiteModels.calc_alpha(v_app, vec_z)
    @test alpha ≈ -1.091003745821884
end

@testset "test_set_cl_cd       " begin
    alpha = deg2rad(10.0)
    KiteModels.set_cl_cd!(kps, alpha)
    @test kps.param_cl ≈ 0.5740976324353215
    @test kps.param_cd ≈ 0.12533689264639403
end

@testset "test_calc_set_cl_cd  " begin
    kps.kcu.depower = 0.236
    kps.kcu.set_depower = kps.kcu.depower
    KiteModels.set_depower_steering!(kps, get_depower(kps.kcu), get_steering(kps.kcu))
    v_app = KVec3(10,2,3)
    vec_c = KVec3(3,2,0)
    KiteModels.calc_set_cl_cd!(kps, vec_c, v_app)
    @test kps.param_cl ≈ -0.09238222805380593
    @test kps.param_cd ≈ 0.8117345278100984
end

@testset "test_clear           " begin
    kps.t_0 = 10.0
    KiteModels.clear!(kps)
    @test kps.t_0 == 0.0
    @test kps.v_reel_out == 0.0
    @test kps.last_v_reel_out == 0.0
end

# Inputs:
# State vector state_y   = pos1, pos2, ..., posn, vel1, vel2, ..., veln
# Derivative   der_yd    = vel1, vel2, ..., veln, acc1, acc2, ..., accn
# Output:
# Residual     res = res1, res2 = pos1,  ..., vel1, ...
@testset "test_residual!       " begin
    res1 = zeros(SVector{SEGMENTS, KVec3})
    res2 = deepcopy(res1)
    res = reduce(vcat, vcat(res1, res2))
    if USE_WINCH
        res = vcat(res, zeros(2))
    end
    X = zeros(SimFloat, 2*kps.set.segments)
    y0, yd0 = KiteModels.init(kps, X)
    # println(y0)
    # println(yd0)
    p = kps
    t = 0.0
    clear!(kps)
    residual!(res, yd0, y0, p, t)
    res1 = res[1:3*SEGMENTS]
    res2 = res[3*SEGMENTS+1:end]
    @test res1 == zeros(3*(SEGMENTS))
    # TODO: add test for res2
    # println(res2)
end

@testset "test_set_v_reel_out  " begin
    v_reel_out = 1.1
    t_0 = 5.5
    KiteModels.set_v_reel_out!(kps, v_reel_out, t_0)
    @test_broken kps.v_reel_out ≈ 1.1
    @test kps.t_0 ≈ 5.5
    clear!(kps)
end

@testset "test_set_depower_steering" begin
    depower  = 0.25
    steering = 0.1
    KiteModels.set_depower_steering!(kps, depower, steering)
    @test kps.depower == depower
    @test kps.steering ≈ 0.09034121603653548
end

@testset "test_init            " begin
    my_state = deepcopy(kps)
    y0, yd0 = KiteModels.init(my_state, zeros(SimFloat, 2*SEGMENTS), delta=1e-6)
    if ! USE_WINCH
        @test length(y0)  == (SEGMENTS) * 6
        @test length(yd0) == (SEGMENTS) * 6
        @test sum(y0)  ≈ 717.163369868302
    else
        @test length(y0)  == (SEGMENTS) * 6 + 2
        @test length(yd0) == (SEGMENTS) * 6 + 2
        @test sum(y0[1:end-2])  ≈ 717.163369868302
    end
    @test sum(yd0) ≈ 3.6e-5
    @test isapprox(my_state.param_cl, 0.574103590856, atol=1e-4)
    @test isapprox(my_state.param_cd, 0.125342896308, atol=1e-4)
end

function test_initial_condition(params::Vector)
    my_state = kps
    y0, yd0 = KiteModels.init(my_state, params)
    residual!(res3, yd0, y0, kps, 0.0)
    return norm(res3) # z component of force on all particles but the first
end

res = nothing
x= nothing
z= nothing
@testset "test_initial_residual" begin
    # global res, x, z
    init_392()
    initial_x =  [-1.52505,  -3.67761,  -5.51761,  -6.08916,  -4.41371,  0.902124,  0.366393,  0.909132,  1.27537,  1.1538,  0.300657,  -1.51768]
    res = test_initial_condition(initial_x)

    my_state = kps
    kps.set.l_tether = 392.0
    kps.set.elevation = 70.0
    kps.set.area = 10.0
    kps.set.v_wind = 9.1
    kps.set.mass = 6.2
    KiteModels.clear!(my_state)
    alpha = deg2rad(10.0)
    KiteModels.set_cl_cd!(kps, alpha)

    # println("state.param_cl: $(my_state.param_cl), state.param_cd: $(my_state.param_cd)")
    # println("res2: "); display(my_state.res2)
    # println("pos: "); display(my_state.pos)
    x = Float64[] 
    z = Float64[]
    for i in 1:length(my_state.pos)
        push!(x, my_state.pos[i][1])
        push!(z, my_state.pos[i][3])
    end  
    # println(norm(res))

    @test my_state.segment_length ≈ 65.33333333333333
    @test my_state.c_spring ≈ 9407.142857142859
    @test my_state.damping  ≈  14.479591836734695
    # @test isapprox(my_state.param_cl, 1.0641931441572074, atol=1e-4)
    # @test isapprox(my_state.param_cd, 0.22825898470541978, atol=1e-4)
    @test sum(my_state.res1) ≈ [0.0, 0.0, 0.0]
    @test my_state.res2[1]   ≈ [0.0, 0.0, 0.0]

    # @test isapprox(my_state.res2[2], [8.83559075e+00, -4.72588546e-07, -5.10109289e+00], rtol=3e-2)
    # @test isapprox(my_state.res2[3], [8.81318565e+00, -4.68864292e-07, -5.08829453e+00], rtol=1e-3)
    # @test isapprox(my_state.res2[7], [1.49735632e+01,  2.71870215e-06,  4.51115984e+01], rtol=1e-3)
    # println("res2: "); display(my_state.res2)
    # println("lift force: $(norm(my_state.lift_force)) N")
end

@testset "test_getters" begin
    x, y, z = kite_ref_frame(kps)
    @test all(x .≈ [-0.9070010101306292, 0.0, 0.4211284455151642])
    @test all(y .≈ [0.0, 1.0, 0.0])
    @test all(z .≈ [-0.4211284455151642, -0.0, -0.9070010101306292])
    @test all(orient_euler(kps) .≈ [1.5707963267948966, -0.4346891114736793, 1.5707963267948966])
    @test all(pos_kite(kps) .≈ [134.97402018366216, 0.0, 366.8418273480761])
    @test calc_elevation(kps) .≈ 1.2182337959242815 # 69.8 deg
    @test calc_azimuth(kps) ≈ 0
    @test calc_heading(kps) ≈ 0
    calc_course(kps) # the course for vel_kite=zero is undefined, so we cannot test it
end

@testset "test_find_steady_state" begin
   KiteModels.set_depower_steering!(kps, 0.25, 0.0)
   res1, res2 = find_steady_state!(kps; delta=1e-6, prn=false) 
   @test norm(res2) < 1e-5                            # velocity and acceleration must be near zero
   pre_tension = KiteModels.calc_pre_tension(kps)
   @test pre_tension > 1.0001
   @test pre_tension < 1.01
   @test unstretched_length(kps) ≈ 392.0              # initial, unstreched tether lenght
   @test tether_length(kps) ≈ 392.1861381318156 # real, streched tether length
   @test winch_force(kps) ≈ 276.25751212763817        # initial force at the winch [N]
   lift, drag = lift_drag(kps)
   @test lift ≈ 443.63277537186394                    # initial lift force of the kite [N]
   @test drag ≈ 94.25218065939362                     # initial drag force of the kite [N]
   @test lift_over_drag(kps) ≈ 4.706870146326417      # initial lift-over-drag
   @test norm(v_wind_kite(kps)) ≈ 9.107670173739065   # initial wind speed at the height of the kite [m/s]
end

function run_benchmarks()
    println("\ncalc_rho:")
    show(@benchmark calc_rho(height) setup=(height=1.0 + rand() * 200.0))
    println("\ncalc_wind_factor:")
    show(@benchmark calc_wind_factor(state.am, height) setup=(height=rand() * 200.0))
    println("\ncalc_cl:")
    show(@benchmark calc_cl(α) setup=(α=(rand()-0.5) * 360.0))
    println("\ncalc_drag:")
    show(@benchmark KiteModels.calc_drag(state, v_segment, unit_vector, rho, v_app_perp, 
                            area) setup=(v_segment = KVec3(1.0, 2, 3);
                            unit_vector = KVec3(2.0, 3.0, 4.0); 
                            rho = calc_rho(10.0f0); state.last_tether_drag = KVec3(0.0, 0.0, 0.0); 
                            v_app_perp =  KVec3(0, -3.0, -4.0); area=se().area))
    println("\ncalc_aero_forces:")
    show(@benchmark KPS3.calc_aero_forces(state, pos_kite, v_kite, rho, rel_steering) setup=(state.v_apparent .= KVec3(35.1,
                                        52.2, 69.3); pos_kite = KVec3(30.0, 5.0, 100.0);  
                                        v_kite = KVec3(3.0, 5.0, 2.0);  
                                        rho = SimFloat(calc_rho(10.0));  rel_steering = 0.1))
    println("\ncalc_res:")
    show(@benchmark KPS3.calc_res(state, pos1, pos2, vel1, vel2, mass, veld, result, i) setup=(i = 1; 
                            pos1 = KVec3(30.0, 5.0, 100.0); pos2 = KVec3(30.0+10, 5.0+11, 100.0+20); 
                            vel1 = KVec3(3.0, 5.0, 2.0); vel2 = KVec3(3.0+0.1, 5.0+0.2, 2.0+0.3); 
                            mass = 9.0; veld = KVec3(0.1, 0.3, 0.4); result = KVec3(0, 0, 0)))
    println("\ncalc_loop")
    show(@benchmark KPS3.loop(state, pos, vel, posd, veld, res1, res2) setup=(pos = zeros(SVector{SEGMENTS+1, KVec3}); 
                            vel  = zeros(SVector{SEGMENTS+1, KVec3}); posd  = zeros(SVector{SEGMENTS+1, KVec3}); 
                            veld  = zeros(SVector{SEGMENTS+1, KVec3}); res1  = zeros(SVector{SEGMENTS+1, KVec3}); 
                            res2  = zeros(SVector{SEGMENTS+1, KVec3}) ))
    println("\nset_cl_cd")
    show(@benchmark KPS3.set_cl_cd!(state, alpha) setup= (alpha = 10.0))
    println("\ncalc_alpha")
    show(@benchmark KPS3.calc_alpha(v_app, vec_z) setup=(v_app = KVec3(10,2,3); vec_z = normalize(KVec3(3,2,0))))
    println("\ncalc_set_cl_cd")
    show(@benchmark KPS3.calc_set_cl_cd!(state, vec_c, v_app) setup=(v_app = KVec3(10,2,3); vec_c = KVec3(3,2,0)))
    println("\nresidual!")
    show(@benchmark residual!(res, yd, y, p, t) setup = (res1 = zeros(SVector{SEGMENTS, KVec3}); res2 = deepcopy(res1); 
                                                            res = reduce(vcat, vcat(res1, res2)); pos = deepcopy(res1);
                                                            pos[1] .= [1.0,2,3]; vel = deepcopy(res1); y = reduce(vcat, vcat(pos, vel));
                                                            der_pos = deepcopy(res1); der_vel = deepcopy(res1); yd = reduce(vcat, vcat(der_pos, der_vel));
                                                            p = SciMLBase.NullParameters(); t = 0.0))
    println("\ntest_initial_condition")
    show(@benchmark res=test_initial_condition(initial_x) setup = (initial_x = (zeros(12))))
    println()
end

end
nothing
