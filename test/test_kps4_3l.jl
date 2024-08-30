using Test, BenchmarkTools, StaticArrays, LinearAlgebra, KiteUtils
using KiteModels, KitePodModels

if ! @isdefined kcu
    const kcu = KCU(se("system_3l.yaml"))
end
kcu.set.winch_model = "AsyncMachine"
kps4_3l::KPS4_3L = KPS4_3L(kcu)

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
    kps4_3l.set.c_s = 2.59
    kps4_3l.set.mass = 0.9
    kps4_3l.set.drum_radius = 0.11
    kps4_3l.set.gear_ratio = 1.0
    kps4_3l.set.inertia_total = 0.104
    kps4_3l.set.profile_law = 3
    kps4_3l.set.sim_settings = "3l_settings.yaml"
    kps4_3l.set.sim_time = 100.0
    kps4_3l.set.abs_tol = 0.0006
    kps4_3l.set.rel_tol = 0.001
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

global integrator
@testset "test_init         " begin
    set_defaults()
    [kps4_3l.pos[i] .= 0 for i in 1:kps4_3l.num_A]

    initial_pos = [
        [-6.1570560327528935e-25 0.0 1.4806772640368532e-22]
        [-2.566743228130072e-24 0.0 -8.086874686735512e-22]
        [4.083892424602851e-22 2.6087432417697598e-23 -5.451085531899676e-21]
        [3.1373804842299133 0.10233230375465678 8.419239953782231]
        [3.137380484229934 -0.10233230375465678 8.419239953782226]
        [2.79081554462356 -2.1852763751589947e-15 7.865164435123338]
        [6.121320828420606 0.20522889118254548 16.89407047913194]
        [6.121320828420612 -0.20522889118254545 16.894070479131944]
        [5.5494165098525565 -8.659299309026002e-15 15.74168566548707]
        [9.08993398697246 0.3081839333650462 25.374281996605426]
        [9.089933986972474 -0.3081839333650463 25.37428199660543]
        [8.299215524060237 -2.368787681503409e-14 23.621283904762933]
        [12.054504391686708 0.41115623825471975 33.85590809205806]
        [12.054504391686724 -0.4111562382547196 33.85590809205806]
        [11.080504250861264 -4.471375328708978e-14 31.489821495039273]
        [15.0171682525174 0.5141319234876086 42.33820094396125]
        [15.017168252517587 -0.514131923487607 42.33820094396119]
        [13.833549993696346 -5.975698338719354e-14 39.368283904135595]
        [17.708219238742384 0.6180997970322778 50.91052388103372]
        [17.70821923874239 -0.6180997970322514 50.91052388103372]
        [16.55629897859964 -6.339537331068523e-14 47.25726620765363]
        [17.73106532237422 0.6180997970322794 50.982979096227474]
        [17.731065322374235 -0.6180997970322498 50.9829790962275]
        [18.49489498763397 1.8603224630271326e-14 50.73718195698237]
    ]

    # initial init
    kps4_3l.set.mass = 0.9
    kps4_3l.set.l_tether = 50.0
    time1 = @elapsed global integrator = KiteModels.init_sim!(kps4_3l; stiffness_factor=1.0, prn=true, torque_control=false)
    pos1 = deepcopy(kps4_3l.pos)
    for i in eachindex(pos1)
        # println(pos1[i]')
        @test isapprox(pos1[i], initial_pos[i,:], atol=1e-5)
    end
    initial_pos = pos1

    # init after changing settings
    kps4_3l.set.mass = 1.0
    kps4_3l.set.l_tether = 51.0
    time2 = @elapsed global integrator = KiteModels.init_sim!(kps4_3l; stiffness_factor=1.0, prn=true, torque_control=false)
    pos2 = deepcopy(kps4_3l.pos)
    @test isapprox(kps4_3l.tether_lengths[3], 51.0, atol=0.1)
    for i in 4:kps4_3l.num_A
        @test all(pos2[i] .!= initial_pos[i])
    end

    # init after changing settings back
    kps4_3l.set.mass = 0.9
    kps4_3l.set.l_tether = 50.0
    time3 = @elapsed global integrator = KiteModels.init_sim!(kps4_3l; stiffness_factor=1.0, prn=true, torque_control=false)
    pos3 = deepcopy(kps4_3l.pos)
    for i in eachindex(pos1)
        @test all(pos3[i] .== initial_pos[i])
    end
    @test isapprox(time2, time3, atol=1.0)

    # init after changing only initial conditions
    kps4_3l.set.elevation = 80.0
    time4 = @elapsed global integrator = KiteModels.init_sim!(kps4_3l; stiffness_factor=1.0, prn=true, torque_control=false)
    pos4 = deepcopy(kps4_3l.pos)
    @test isapprox(rad2deg(calc_elevation(kps4_3l)), 80.0, atol=2.0)
    @test time3 / time4 > 5.0
    for i in 4:kps4_3l.num_A
        @test all(pos4[i] .!= initial_pos[i])
    end

    # init after just stepping
    KiteModels.next_step!(kps4_3l, integrator)
    time5 = @elapsed global integrator = KiteModels.init_sim!(kps4_3l; stiffness_factor=1.0, prn=true, torque_control=false)
    pos5 = deepcopy(kps4_3l.pos)
    @test time3 / time5 > 100
    for i in eachindex(pos1)
        @test all(pos5[i] .== pos4[i])
    end
end

@testset "test_step         " begin
    kps4_3l.set.elevation = 70.8
    time4 = @elapsed global integrator = KiteModels.init_sim!(kps4_3l; stiffness_factor=1.0, prn=true, torque_control=false)

    # KiteModels.next_step!(kps4_3l, integrator)
    pos2 = [
        [-6.1570560327528935e-25 0.0 1.4806772640368532e-22]
        [-2.566743228130072e-24 0.0 -8.086874686735512e-22]
        [4.083892424602851e-22 2.6087432417697598e-23 -5.451085531899676e-21]
        [3.1373804842299133 0.10233230375465678 8.419239953782231]
        [3.137380484229934 -0.10233230375465678 8.419239953782226]
        [2.79081554462356 -2.1852763751589947e-15 7.865164435123338]
        [6.121320828420606 0.20522889118254548 16.89407047913194]
        [6.121320828420612 -0.20522889118254545 16.894070479131944]
        [5.5494165098525565 -8.659299309026002e-15 15.74168566548707]
        [9.08993398697246 0.3081839333650462 25.374281996605426]
        [9.089933986972474 -0.3081839333650463 25.37428199660543]
        [8.299215524060237 -2.368787681503409e-14 23.621283904762933]
        [12.054504391686708 0.41115623825471975 33.85590809205806]
        [12.054504391686724 -0.4111562382547196 33.85590809205806]
        [11.080504250861264 -4.471375328708978e-14 31.489821495039273]
        [15.0171682525174 0.5141319234876086 42.33820094396125]
        [15.017168252517587 -0.514131923487607 42.33820094396119]
        [13.833549993696346 -5.975698338719354e-14 39.368283904135595]
        [17.708219238742384 0.6180997970322778 50.91052388103372]
        [17.70821923874239 -0.6180997970322514 50.91052388103372]
        [16.55629897859964 -6.339537331068523e-14 47.25726620765363]
        [17.73106532237422 0.6180997970322794 50.982979096227474]
        [17.731065322374235 -0.6180997970322498 50.9829790962275]
        [18.49489498763397 1.8603224630271326e-14 50.73718195698237]
    ]
    for i in eachindex(kps4_3l.pos)
        # println(kps4_3l.pos[i]')
        @test isapprox(pos2[i,:], kps4_3l.pos[i], atol=1e-6)
    end
    # println(kps4_3l.L_C)
    @test all(kps4_3l.L_C .≈ [22.65112017603021, 162.37893030596314, 466.36470393015884])
end

function simulate(integrator, steps)
    for i in 1:steps
        KiteModels.next_step!(kps4_3l, integrator; set_values=[0.0, 0.0, 0.15])
        @show i
    end
    return integrator.iter/steps
end

@testset "test_simulate     " begin
    STEPS = 20
    integrator = KiteModels.init_sim!(kps4_3l; stiffness_factor=1.0, prn=true, torque_control=false)
    # println("\nStarting simulation...")
    av_steps = simulate(integrator, STEPS)
    if Sys.isapple()
        println("isapple $av_steps")
        @test isapprox(av_steps, 10.85, atol=1.0)
    else
        println("not apple $av_steps")
        @test isapprox(av_steps, 10.85, atol=1.0)
    end
  
    @test -10.0 < kps4_3l.L_C[1] < 10.0
    @test 150.0 < kps4_3l.L_C[2] < 200.0
    @test 400.0 < kps4_3l.L_C[3] < 600.0
    @test isapprox(kps4_3l.reel_out_speeds, [0.15824099721234128, 0.15824112269822727, 0.20901760480448708], atol=1e-6)
    @test isapprox(kps4_3l.L_C[2], -kps4_3l.L_D[2], atol=1e-2)
    @show kps4_3l.L_C
    @show kps4_3l.L_D
    @show kps4_3l.reel_out_speeds
    
    # TODO Add testcase with varying reelout speed 
end

# TODO: add testset for sysstate

end
nothing