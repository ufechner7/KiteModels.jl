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
    KiteModels.clear!(kps4_3l)
    # kps4_3l.set.
end

set_defaults()

pos1 = nothing
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
    global pos1 = deepcopy(kps4_3l.pos)
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
    update_pos!(kps4_3l, integrator)
    @test all(kps4_3l.pos[4:kps4_3l.num_A] .!= pos1[4:kps4_3l.num_A])
    @test integrator.last_stepfail == false
end

@testset "test_init         " begin
    
end

@testset "test_step         " begin
    
end

@testset "test_reset        " begin
    reset_sim!
end

@testset "test_simulate     " begin
    # STEPS = 50

    # kps4_3l.set.solver = "DFBDF"
    # # println("finding steady state")
    # init_50()
    # integrator = KiteModels.init_sim!(kps4_3l; stiffness_factor=0.035, prn=false)
    # # println("\nStarting simulation...")
    # simulate(integrator, STEPS)
    # av_steps = simulate(integrator, STEPS-10)
    # if Sys.isapple()
    #     println("isapple $av_steps")
    #     @test isapprox(av_steps, 835.25, rtol=0.6)
    # else
    #     println("not apple $av_steps")
    #     @test isapprox(av_steps, 835.25, rtol=0.6)
    # end
  
    # lift, drag = KiteModels.lift_drag(kps4_3l)
    # # println(lift, " ", drag) # 703.7699568972286 161.44746368100536
    # @test isapprox(lift, 404.2596735903995, rtol=0.05)
    # sys_state = SysState(kps4_3l)
    # update_sys_state!(sys_state, kps4_3l)
    # # TODO Add testcase with varying reelout speed 
end
# println(kps4_3l.set)

end
nothing