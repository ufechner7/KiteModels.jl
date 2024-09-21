using Test, BenchmarkTools, StaticArrays, LinearAlgebra, KiteUtils
using KiteModels, KitePodModels

kcu_3l::KCU = KCU(se("system_3l.yaml"))
kcu_3l.set.winch_model = "AsyncMachine"
kps4_3l::KPS4_3L = KPS4_3L(kcu_3l)

pos, vel = nothing, nothing

@testset verbose = true "KPS4_3L tests...." begin

function set_defaults()
    kps4_3l.set = update_settings()
    KiteModels.clear!(kps4_3l)
end

set_defaults()

@testset "test_init         " begin
    set_defaults()
    [kps4_3l.pos[i] .= 0 for i in 1:kps4_3l.num_A]

    initial_pos = [
        [4.267845479738403e-15 1.1257108306179438e-18 -2.5589696375556316e-15]
        [-4.726038177476628e-15 -7.146085596181844e-18 -6.539688409891346e-16]
        [-1.7852893341457758e-16 1.3427794793615233e-16 1.2233052837128663e-16]
        [7.8512708081282065 0.10807103464352219 6.8857661642582135]
        [7.85127080633873 -0.1080710349791022 6.885766166289688]
        [0.31547803473993047 -2.551069830707611e-10 8.361548956869619]
        [12.550117120853807 0.24785043130714346 16.21150568397882]
        [12.550117121787169 -0.24785043177550484 16.211505684636375]
        [0.6195651853664321 -5.144368422678874e-10 16.723528875167748]
        [14.587760449784694 0.3970698904122888 26.453286171590342]
        [14.587760453224249 -0.39706989124115555 26.453286171745447]
        [0.909834125470212 -7.828208252535404e-10 25.086008808615002]
        [13.011486003441334 0.5403517370214007 36.7762433967783]
        [13.011486008461219 -0.540351738478159 36.77624339716316]
        [1.1846892780717257 -1.0662883772808636e-9 33.449018469129555]
        [9.423156013556122 0.6632379993242858 46.58327549727098]
        [9.42315601222154 -0.663238001756424 46.58327549531171]
        [1.4428086869274315 -1.3729954508277103e-9 41.812570299360196]
        [1.9978661022246724 0.7611102069744416 53.92686686299797]
        [1.9978661026765974 -0.7611102109569202 53.92686686282812]
        [1.6830056646333187 -1.7145184505046268e-9 50.17666491922819]
        [1.7273030489131658 0.7979882125246655 54.01496139682041]
        [1.7273030493930808 -0.7979882166941015 54.01496139665826]
        [2.3005532222358047 -1.910461039103992e-9 54.00300139171894]
    ]

    # initial init
    kps4_3l.set.mass = 0.9
    kps4_3l.set.l_tether = 50.0
    KiteModels.init_sim!(kps4_3l; prn=true, torque_control=false)
    pos1 = deepcopy(kps4_3l.pos)
    for i in eachindex(pos1)
        # println(pos1[i]')
        @test isapprox(pos1[i], initial_pos[i,:], atol=1e-2)
    end
    initial_pos = pos1

    # init after changing settings
    kps4_3l.set.mass = 1.0
    kps4_3l.set.l_tether = 51.0
    KiteModels.init_sim!(kps4_3l; prn=true, torque_control=false)
    pos2 = deepcopy(kps4_3l.pos)
    @test isapprox(kps4_3l.tether_lengths[3], 51.0, atol=0.1)
    for i in 4:kps4_3l.num_A
        @test !isapprox(pos2[i], initial_pos[i], atol=1e-2)
    end

    # init after changing settings back
    kps4_3l.set.mass = 0.9
    kps4_3l.set.l_tether = 50.0
    KiteModels.init_sim!(kps4_3l; prn=true, torque_control=false)
    pos3 = deepcopy(kps4_3l.pos)
    for i in eachindex(pos1)
        @test all(pos3[i] .== initial_pos[i])
    end

    # # init after changing only initial conditions
    # kps4_3l.set.elevation = 80.0
    # KiteModels.init_sim!(kps4_3l; prn=true, torque_control=false)
    # pos4 = deepcopy(kps4_3l.pos)
    # @test isapprox(rad2deg(calc_elevation(kps4_3l)), 80.0, atol=2.0)
    # for i in 4:kps4_3l.num_A
    #     @test all(pos4[i] .!= initial_pos[i])
    # end

    # # init after just stepping
    # KiteModels.next_step!(kps4_3l)
    # KiteModels.init_sim!(kps4_3l; prn=true, torque_control=false)
    # pos5 = deepcopy(kps4_3l.pos)
    # for i in eachindex(pos1)
    #     @test all(pos5[i] .== pos4[i])
    # end

    # TODO: add tests for torque controlled
end

# @testset "test_step         " begin
#     kps4_3l.set.elevation = 70.8
#     KiteModels.init_sim!(kps4_3l; prn=true, torque_control=false)

#     # KiteModels.next_step!(kps4_3l)
#     pos2 = [
#         [-6.1570560327528935e-25 0.0 1.4806772640368532e-22]
#         [-2.566743228130072e-24 0.0 -8.086874686735512e-22]
#         [4.083892424602851e-22 2.6087432417697598e-23 -5.451085531899676e-21]
#         [3.1373804842299133 0.10233230375465678 8.419239953782231]
#         [3.137380484229934 -0.10233230375465678 8.419239953782226]
#         [2.79081554462356 -2.1852763751589947e-15 7.865164435123338]
#         [6.121320828420606 0.20522889118254548 16.89407047913194]
#         [6.121320828420612 -0.20522889118254545 16.894070479131944]
#         [5.5494165098525565 -8.659299309026002e-15 15.74168566548707]
#         [9.08993398697246 0.3081839333650462 25.374281996605426]
#         [9.089933986972474 -0.3081839333650463 25.37428199660543]
#         [8.299215524060237 -2.368787681503409e-14 23.621283904762933]
#         [12.054504391686708 0.41115623825471975 33.85590809205806]
#         [12.054504391686724 -0.4111562382547196 33.85590809205806]
#         [11.080504250861264 -4.471375328708978e-14 31.489821495039273]
#         [15.0171682525174 0.5141319234876086 42.33820094396125]
#         [15.017168252517587 -0.514131923487607 42.33820094396119]
#         [13.833549993696346 -5.975698338719354e-14 39.368283904135595]
#         [17.708219238742384 0.6180997970322778 50.91052388103372]
#         [17.70821923874239 -0.6180997970322514 50.91052388103372]
#         [16.55629897859964 -6.339537331068523e-14 47.25726620765363]
#         [17.73106532237422 0.6180997970322794 50.982979096227474]
#         [17.731065322374235 -0.6180997970322498 50.9829790962275]
#         [18.49489498763397 1.8603224630271326e-14 50.73718195698237]
#     ]
#     for i in eachindex(kps4_3l.pos)
#         # println(kps4_3l.pos[i]')
#         @test isapprox(pos2[i,:], kps4_3l.pos[i], atol=1e-6)
#     end
#     # println(kps4_3l.L_C)
#     @test all(kps4_3l.L_C .≈ [22.65112017603021, 162.37893030596314, 466.36470393015884])
# end

# function simulate(integrator, steps)
#     for i in 1:steps
#         KiteModels.next_step!(kps4_3l; set_values=[0.0, 0.0, 0.15])
#     end
#     return integrator.iter/steps
# end

# @testset "test_simulate     " begin
#     STEPS = 20
#     integrator = KiteModels.init_sim!(kps4_3l; prn=true, torque_control=false)
#     # println("\nStarting simulation...")
#     av_steps = simulate(integrator, STEPS)
#     if Sys.isapple()
#         println("isapple $av_steps")
#         @test isapprox(av_steps, 10.85, atol=1.0)
#     else
#         println("not apple $av_steps")
#         @test isapprox(av_steps, 10.85, atol=1.0)
#     end
  
#     @test -10.0 < kps4_3l.L_C[1] < 10.0
#     @test 150.0 < kps4_3l.L_C[2] < 200.0
#     @test 400.0 < kps4_3l.L_C[3] < 600.0
#     @test isapprox(kps4_3l.reel_out_speeds, [0.15824099721234128, 0.15824112269822727, 0.20901760480448708], atol=1e-6)
#     @test isapprox(kps4_3l.L_C[2], -kps4_3l.L_D[2], atol=1e-2)
    
#     # TODO Add testcase with varying reelout speed 
# end

# TODO: add testset for sysstate

end
nothing