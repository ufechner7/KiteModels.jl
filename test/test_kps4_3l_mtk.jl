using Test, BenchmarkTools, StaticArrays, LinearAlgebra, KiteUtils, ModelingToolkit, SymbolicIndexingInterface, OrdinaryDiffEq
using KiteModels, KitePodModels

if ! @isdefined kcu
    const kcu = KCU(se("system_3l.yaml"))
end
if ! @isdefined kps4_3l
    kcu.set.winch_model = "AsyncMachine"
    const kps4_3l = KPS4_3L(kcu)
end

pos, vel = nothing, nothing

@testset verbose = true "KPS4_3L_MTK tests...." begin

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
    kps4_3l.set.c_s = 2.59
    kps4_3l.set.mass = 0.9
    kps4_3l.set.drum_radius = 0.11
    kps4_3l.set.gear_ratio = 1.0
    kps4_3l.set.inertia_total = 0.104
    kps4_3l.set.profile_law = 3
    kps4_3l.set.sim_settings = "3l_settings.yaml"
    kps4_3l.set.sim_time = 100.0
    kps4_3l.set.abs_tol = 0.006
    kps4_3l.set.rel_tol = 0.01
    kps4_3l.set.max_iter = 10000
    kps4_3l.set.physical_model = "KPS4_3L"
    kps4_3l.set.version = 2
    kps4_3l.set.cl_list = [0.0, 0.5, 0.0, 0.08, 0.125, 0.15, 0.0, 1.0, 1.0, 0.0, -0.5, 0.0]
    kps4_3l.set.cd_list = [0.5, 0.5, 0.5,  1.0,   0.2,  0.1, 0.2, 1.0, 0.5, 0.5,  0.5]
    kps4_3l.set.d_tether = 1.0
    kps4_3l.set.v_wind_ref = [15.51, 0.0]
    kps4_3l.set.v_wind = 15.51
    kps4_3l.set.rel_turbs = [0.0, 0.0, 0.0]
    KiteModels.clear!(kps4_3l)
end

set_defaults()

@testset "test_model!       " begin
    kps4_3l.stiffness_factor = 0.04
    res =  zeros(MVector{6*(kps4_3l.num_A-5)+4+6, SimFloat})
    kps4_3l.mtk = true
    y0, _ = KiteModels.find_steady_state!(kps4_3l; stiffness_factor=0.1, prn=true)
    y0_ = [4.490372205878421, 0.10298757652176786, 7.780950480262341, 4.490372205878421, -0.10298757652176786, 7.780950480262341, 4.166666666666668, 1.0e-6, 7.216878364870322, 8.980744411756842, 0.20597515304353572, 15.561900960524682, 8.980744411756842, -0.20597515304353572, 15.561900960524682, 8.333333333333336, 1.0e-6, 14.433756729740644, 13.471116617635262, 0.3089627295653036, 23.342851440787022, 13.471116617635262, -0.3089627295653036, 23.342851440787022, 12.500000000000004, 1.0e-6, 21.650635094610966, 17.961488823513683, 0.41195030608707145, 31.123801921049363, 17.961488823513683, -0.41195030608707145, 31.123801921049363, 16.66666666666667, 1.0e-6, 28.867513459481287, 22.451861029392106, 0.5149378826088393, 38.904752401311704, 22.451861029392106, -0.5149378826088393, 38.904752401311704, 20.83333333333334, 1.0e-6, 36.08439182435161, 25.000000000000007, 1.0e-6, 43.30127018922193, 26.942233235270525, 0.6179254591306071, 46.685702881574045, 26.942233235270525, -0.6179254591306071, 46.685702881574045, 27.646090205747516, 0.013239546805801463, 46.27933087019816, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 0.0, 0.0, 0.0, 0.0, 53.90566405419164, 53.90566405419164, 50.0, 0.0, 0.0, 0.0]
    # println(y0)
    # println(yd0)
    # @test all(y0 .≈ y0_)

    time = 0.0
    dt = 1/kps4_3l.set.sample_freq
    tspan   = (0.0, dt) 
    abstol  = kps4_3l.set.abs_tol 
    @test kps4_3l.num_A == length(y0)
    simple_sys, sys = model!(kps4_3l, y0; torque_control=false)
    @test length(equations(sys)) > length(equations(simple_sys))
    @test length(unknowns(simple_sys)) == length(equations(simple_sys))
    @test length(equations(simple_sys)) == 142
    kps4_3l.prob = ODEProblem(simple_sys, nothing, tspan)
    integrator = OrdinaryDiffEq.init(kps4_3l.prob, KenCarp4(autodiff=false); dt, abstol=kps4_3l.set.abs_tol, reltol=kps4_3l.set.rel_tol, save_on=false)
    @test length(integrator.u) == 142
    kps4_3l.set_values_idx = parameter_index(integrator.f, :set_values)
    kps4_3l.v_wind_gnd_idx = parameter_index(integrator.f, :v_wind_gnd)
    kps4_3l.v_wind_idx = parameter_index(integrator.f, :v_wind)
    kps4_3l.stiffness_factor_idx = parameter_index(integrator.f, :stiffness_factor)
    kps4_3l.get_pos = getu(integrator.sol, simple_sys.pos[:,:])
    kps4_3l.get_steering_pos = getu(integrator.sol, simple_sys.steering_pos)
    kps4_3l.get_line_acc = getu(integrator.sol, simple_sys.acc[:,kps4_3l.num_E-2])
    kps4_3l.get_kite_vel = getu(integrator.sol, simple_sys.vel[:,kps4_3l.num_A])
    kps4_3l.get_winch_forces = getu(integrator.sol, simple_sys.force[:,1:3])
    kps4_3l.get_L_C = getu(integrator.sol, simple_sys.L_C)
    kps4_3l.get_L_D = getu(integrator.sol, simple_sys.L_D)
    kps4_3l.get_D_C = getu(integrator.sol, simple_sys.D_C)
    kps4_3l.get_D_D = getu(integrator.sol, simple_sys.D_D)
    kps4_3l.get_tether_lengths = getu(integrator.sol, simple_sys.tether_length)
    kps4_3l.get_tether_speeds = getu(integrator.sol, simple_sys.tether_speed)
    update_pos!(kps4_3l, integrator)

    for (p, y) in zip(kps4_3l.pos, y0)
        @test all(p .== y)
    end
    @test all(kps4_3l.steering_pos .== 0)
    @test all(kps4_3l.vel_kite .== 0)
    @test isapprox(kps4_3l.L_C, kps4_3l.L_D .* [1,-1,1], atol=1e-5)
    @test isapprox(kps4_3l.D_C, kps4_3l.D_D .* [1,-1,1], atol=1e-5)
    @test all(kps4_3l.L_C .≈ [-1.050236194673068e-16, 112.36591385162515, 327.9801594516224])
    @test all(kps4_3l.D_C .≈ [70.49414103359018, 4.935335969053079, -2.364032694118963])
    # println(kps4_3l.L_C)
    # println(kps4_3l.D_C)


    # test step
    pos1 = deepcopy(kps4_3l.pos)
    kps4_3l.stiffness_factor = 1.0
    kps4_3l.iter = 0
    KiteModels.set_v_wind_ground!(kps4_3l, calc_height(kps4_3l), kps4_3l.set.v_wind, 0.0)
    kps4_3l.set_values .= [0.4, 0.5, 0.6]
    integrator.ps[kps4_3l.set_values_idx] .= kps4_3l.set_values
    integrator.ps[kps4_3l.v_wind_gnd_idx] .= kps4_3l.v_wind_gnd
    integrator.ps[kps4_3l.v_wind_idx] .= kps4_3l.set.v_wind
    integrator.ps[kps4_3l.stiffness_factor_idx] = kps4_3l.stiffness_factor
    @test all(integrator.ps[kps4_3l.set_values_idx] .== kps4_3l.set_values)
    kps4_3l.t_0 = integrator.t
    OrdinaryDiffEq.step!(integrator, dt, true)
    @test all(kps4_3l.pos[4:kps4_3l.num_A] .== pos1[4:kps4_3l.num_A])
    update_pos!(kps4_3l, integrator)
    @test all(kps4_3l.pos[4:kps4_3l.num_A] .!= pos1[4:kps4_3l.num_A])
    @test integrator.last_stepfail == false
    @test kps4_3l.mtk == true
end

global integrator, pos1
@testset "test_init         " begin
    set_defaults()
    [kps4_3l.pos[i] .= 0 for i in 1:kps4_3l.num_A]
    global integrator = KiteModels.init_sim!(kps4_3l; stiffness_factor=0.1, prn=false, mtk=true, torque_control=false)
    @test kps4_3l.mtk == true
    global pos1 = [
        [0.0 0.0 0.0]
        [0.0 0.0 0.0]
        [0.0 0.0 0.0]
        [2.826555134916249 0.10162113192687891 8.622005105305957]
        [2.826555134916249 -0.10162113192687891 8.622005105305957]
        [2.1298311552238363 0.0 8.131298778702757]
        [5.456422706641912 0.2038737369459307 17.306023661511617]
        [5.456422706641912 -0.2038737369459307 17.306023661511617]
        [4.241908474517871 0.0 16.26723340006626]
        [7.850932935818464 0.3068255041595917 26.05786267666516]
        [7.850932935818464 -0.3068255041595917 26.05786267666516]
        [6.333100132564329 0.0 24.408567603645043]
        [9.98601991466058 0.4104649244760199 34.87657889907129]
        [9.98601991466058 -0.4104649244760199 34.87657889907129]
        [8.401700377570428 0.0 32.55567756458501]
        [11.842878698197364 0.5147433105579666 43.758044143291855]
        [11.842878698197364 -0.5147433105579666 43.758044143291855]
        [10.446505443256301 0.0 40.70879859662597]
        [13.405815986564694 0.6195842015236318 52.69591805323245]
        [13.405815986564694 -0.6195842015236318 52.69591805323245]
        [12.466577313446024 0.0 48.868089171558495]
        [13.405815986564694 0.6195842015236318 52.69591805323245]
        [13.405815986564694 -0.6195842015236318 52.69591805323245]
        [14.171760968247852 1.1102230246251565e-16 52.4606698480763]
    ]
    for i in eachindex(kps4_3l.pos)
        # println(kps4_3l.pos[i]')
        @test all(pos1[i,:] .≈ kps4_3l.pos[i])
    end
    sys_state = KiteModels.SysState(kps4_3l)

end

@testset "test_step         " begin
    @test kps4_3l.mtk == true
    KiteModels.next_step!(kps4_3l, integrator)
    pos2 = [
        [-5.796369547658875e-19 3.105564228811931e-21 2.048931634380671e-18]
        [4.1775930164545956e-20 1.2010712923916298e-21 7.298750574041877e-19]
        [-1.0572099515026155e-18 3.3782511637807305e-19 5.274874242660198e-18]
        [2.8298579736741147 0.10161329821195392 8.621175820966815]
        [2.830623570190892 -0.10161108533187954 8.620950345308996]
        [2.1302535909219413 -1.292341650237965e-5 8.132038539148022]
        [5.460582249504757 0.20386598364954833 17.305187161613578]
        [5.461538411996239 -0.20386344965793496 17.304929609827887]
        [4.242016191052998 -4.1371222398968844e-5 16.26890498869265]
        [7.855454213540059 0.3068196645628308 26.05717753387838]
        [7.856489915153269 -0.30681716661732195 26.056923670495657]
        [6.3316955459341155 -0.00010041801729104318 24.411477480427358]
        [9.99084415164687 0.4104610182125188 34.87606887722985]
        [9.991941232468655 -0.41045821413403777 34.8758254139833]
        [8.397847155693308 -0.00018824714259965992 32.56005879086417]
        [11.846400321894027 0.5147447186417542 43.75805306941343]
        [11.847471812019844 -0.5147291067890325 43.75784017580861]
        [10.440210383536016 -0.0002653046876820298 40.71464210076245]
        [13.393727106942748 0.619648084377197 52.69888711431324]
        [13.393746698135939 -0.6195315900226024 52.698882288898126]
        [12.457316981603986 -0.0003032160793861237 48.87551692338118]
        [13.394882292641316 0.619648530285251 52.703603694641174]
        [13.394931457320247 -0.6195311326990197 52.703719616670604]
        [14.161416203376493 6.691191585922139e-5 52.47012644562449]
    ]
    for i in eachindex(kps4_3l.pos)
        # println(kps4_3l.pos[i]')
        @test all(pos1[i,:] .!= kps4_3l.pos[i])
        @test isapprox(pos2[i,:], kps4_3l.pos[i], atol=0.005) # TODO: somehow slightly different values when running runtests.jl
    end
    # println("L_C ", kps4_3l.L_C)
    @test isapprox(kps4_3l.L_C, [3.966773080945373, 115.50991115735522, 335.6996775000222], rtol=0.01)
end

@testset "test_reset        " begin
    reset_sim!(kps4_3l)
    reset_time = @elapsed reset_sim!(kps4_3l)
    @test reset_time < 0.01
    for i in eachindex(kps4_3l.pos)
        # println(kps4_3l.pos[i]')
        @test all(pos1[i,:] .== kps4_3l.pos[i])
    end
    @test isapprox(kps4_3l.L_C[1], 0.0, atol=1e-6)
end

function simulate(integrator, steps)
    for i in 1:steps
        KiteModels.next_step!(kps4_3l, integrator; set_values=[0.0, 0.0, 0.35])
    end
    return integrator.iter/steps
end

@testset "test_simulate     " begin
    STEPS = 10
    reset_sim!(kps4_3l)
    # println("\nStarting simulation...")
    simulate(integrator, STEPS)
    av_steps = simulate(integrator, STEPS)
    if Sys.isapple()
        println("isapple $av_steps")
        @test isapprox(av_steps, 11.5, atol=1.0)
    else
        println("not apple $av_steps")
        @test isapprox(av_steps, 11.5, atol=1.0)
    end
  
    println(kps4_3l.L_C)
    @test isapprox(kps4_3l.L_C, [11.298536148915304, 212.2336998177928, 595.3408008967488], atol=0.1)
    
    # @test (normalize(kps4_3l.e_z) - normalize(kps4_3l.L_C))[1] > 0
    # # println(lift, " ", drag) # 703.7699568972286 161.44746368100536
    # @test isapprox(lift, 404.2596735903995, rtol=0.05)
    # sys_state = SysState(kps4_3l)
    # update_sys_state!(sys_state, kps4_3l)
    # # TODO Add testcase with varying reelout speed 
end
# println(kps4_3l.set)

end
nothing