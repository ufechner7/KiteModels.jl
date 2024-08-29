using Test, BenchmarkTools, StaticArrays, LinearAlgebra, KiteUtils, ModelingToolkit, SymbolicIndexingInterface, OrdinaryDiffEq
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

global integrator, initial_pos
@testset "test_init         " begin
    set_defaults()
    [kps4_3l.pos[i] .= 0 for i in 1:kps4_3l.num_A]

    time1 = @elapsed global integrator = KiteModels.init_sim!(kps4_3l; stiffness_factor=1.0, prn=true, torque_control=false)
    pos1 = deepcopy(kps4_3l.pos)
    global initial_pos = pos1

    kps4_3l.set.mass = 0.9 + 1e-9
    kps4_3l.set.l_tether = 50.0 + 1e-9
    time2 = @elapsed global integrator = KiteModels.init_sim!(kps4_3l; stiffness_factor=1.0, prn=true, torque_control=false)
    pos2 = deepcopy(kps4_3l.pos)

    kps4_3l.set.l_tether = 50.0
    time3 = @elapsed global integrator = KiteModels.init_sim!(kps4_3l; stiffness_factor=1.0, prn=true, torque_control=false)
    pos3 = deepcopy(kps4_3l.pos)

    time4 = @elapsed global integrator = KiteModels.init_sim!(kps4_3l; stiffness_factor=1.0, prn=true, torque_control=false)
    pos4 = deepcopy(kps4_3l.pos)
    
    @show time1 time2 time3 time4
    @test time3 < time2/2 && time4 < time2/100
    for i in eachindex(pos1)
        @test isapprox(pos1[i], pos2[i], atol=0.2)
        @test isapprox(pos1[i], pos3[i], atol=0.2)
        @test all(pos3[i] .== pos4[i])
    end
    @show pos3 pos4

    set_defaults()
    integrator = KiteModels.init_sim!(kps4_3l; stiffness_factor=1.0, prn=true, torque_control=false)
    for i in eachindex(kps4_3l.pos)
        println(kps4_3l.pos[i]')
        # @test all(initial_pos[i,:] .≈ kps4_3l.pos[i])
    end
    println("acc ", norm(integrator[kps4_3l.simple_sys.acc]))
    # 1064.7050756272429
    sys_state = KiteModels.SysState(kps4_3l)

end

# @testset "test_step         " begin
#     KiteModels.next_step!(kps4_3l, integrator)
#     pos2 = [
#         [6.194391874676442e-21 -6.603695400302634e-24 -8.888702554685316e-18]
#         [-5.290657911512011e-20 4.0239900676426147e-23 -4.000371890774165e-19]
#         [-1.1914048048487081e-17 -7.349718070579275e-19 -4.009352536206954e-18]
#         [3.4397703933935504 0.10117316797602434 8.304450386921632]
#         [3.4397703931753387 -0.1011731679770315 8.304450387020093]
#         [2.8001918730238726 -3.781549653090172e-11 7.905259572908571]
#         [6.544140772357473 0.20367508223264025 16.740009796269447]
#         [6.544140772335404 -0.20367508223302513 16.740009796303717]
#         [5.62866700899026 -4.9300464869613585e-11 15.800445989330704]
#         [9.540042582094188 0.3065665022849803 25.21470046380188]
#         [9.540042582097662 -0.3065665022853301 25.214700463835058]
#         [8.507496499227601 -6.194824851803334e-11 23.677411705877468]
#         [12.528964235460489 0.4094809546323497 33.69186088445923]
#         [12.528964235497488 -0.4094809546306103 33.6918608844886]
#         [11.339993603079177 -8.273042136528275e-11 31.571151782496848]
#         [15.412646362580034 0.5131843477721921 42.2053911086091]
#         [15.412646362363924 -0.5131843477640395 42.20539110873228]
#         [14.133651994703683 -7.382413369684321e-11 39.47871263242362]
#         [17.864482911376157 0.6193429585590792 50.8531363529247]
#         [17.86448291146066 -0.61934295851386 50.85313635297098]
#         [16.863803598833563 -6.874570775236859e-11 47.408414913865585]
#         [17.955580153794777 0.6193429585673957 51.16672795020686]
#         [17.955580153870212 -0.6193429585055434 51.166727950221926]
#         [18.72161367387903 7.489713534546274e-11 50.93184293581413]
#     ]
#     for i in eachindex(kps4_3l.pos)
#         # println(kps4_3l.pos[i]')
#         @test all(pos2[i,:] .≈ kps4_3l.pos[i])
#     end
#     # println(kps4_3l.L_C)
#     @test all(kps4_3l.L_C .≈ [35.23188440838821, 123.05256310087925, 334.9023917920065])
# end

# function simulate(integrator, steps)
#     for i in 1:steps
#         KiteModels.next_step!(kps4_3l, integrator; set_values=[0.0, 0.0, 0.05])
#     end
#     return integrator.iter/steps
# end

# @testset "test_simulate     " begin
#     STEPS = 30
#     reset_sim!(kps4_3l)
#     # println("\nStarting simulation...")
#     simulate(integrator, STEPS)
#     av_steps = simulate(integrator, STEPS)
#     if Sys.isapple()
#         println("isapple $av_steps")
#         @test isapprox(av_steps, 5.766666666666667, atol=1.0)
#     else
#         println("not apple $av_steps")
#         @test isapprox(av_steps, 5.766666666666667, atol=1.0)
#     end
  
#     # println(kps4_3l.L_C)
#     @test all(kps4_3l.L_C .≈ [-1.577566498953694, 168.94195353362062, 456.44628119609695])
    
#     # @test (normalize(kps4_3l.e_z) - normalize(kps4_3l.L_C))[1] > 0
#     # # println(lift, " ", drag) # 703.7699568972286 161.44746368100536
#     # @test isapprox(lift, 404.2596735903995, rtol=0.05)
#     # sys_state = SysState(kps4_3l)
#     # update_sys_state!(sys_state, kps4_3l)
#     # # TODO Add testcase with varying reelout speed 
# end
# println(kps4_3l.set)

end
nothing