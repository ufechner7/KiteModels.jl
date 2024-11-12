using Test, BenchmarkTools, StaticArrays, LinearAlgebra, KiteUtils
using KiteModels, KitePodModels

old_path = get_data_path()
set_data_path(joinpath(dirname(dirname(pathof(KiteModels))), "data"))
kcu_3l::KCU = KCU(se("system_3l.yaml"))
kcu_3l.set.winch_model = "AsyncMachine"
k3l::KPS4_3L = KPS4_3L(kcu_3l)

pos, vel = nothing, nothing

@testset verbose = true "KPS4_3l tests..." begin

tol::Float32 = 1e-5
prn::Bool = false

function set_defaults()
    k3l.set = update_settings()
    k3l.set.abs_tol = tol
    k3l.set.rel_tol = tol
    KiteModels.clear!(k3l)
end

set_defaults()

global initial_pos
@testset "test_init         " begin
    set_defaults()
    [k3l.pos[i] .= 0 for i in 1:k3l.num_A]

    # initial_pos = [
    #     [1.956416680572584e-18 -3.125513814477895e-21 -5.720339822309706e-19]
    #     [1.91352691093351e-18 4.282925325306576e-20 2.2118754797786398e-18]
    #     [-9.004262075908625e-19 -8.743479042208521e-20 1.043447107448909e-18]
    #     [7.784853189744222 0.10569510168880981 6.960672205291082]
    #     [7.7848531896951885 -0.1056951016892796 6.96067220534591]
    #     [0.6426154451758553 6.949885399751856e-12 8.336738719271432]
    #     [11.586653042156335 0.2506120857714412 16.68648191133713]
    #     [11.586653042075785 -0.2506120857731139 16.686481911404254]
    #     [1.271299088438001 1.3894049639298299e-11 16.674548626145384]
    #     [12.772490415263762 0.4027982287207201 27.061294344744]
    #     [12.772490415114174 -0.4027982287201746 27.061294344819043]
    #     [1.8831522309211308 2.0828036927412173e-11 25.01361953696909]
    #     [12.705233401060338 0.5381159139350681 37.503674710779734]
    #     [12.705233401020347 -0.5381159139285039 37.503674710855556]
    #     [2.4763455426643515 2.7746793461421894e-11 33.354047648677586]
    #     [12.717685940855144 0.6357504779887956 47.94668473725685]
    #     [12.717685940906842 -0.6357504779549034 47.94668473733282]
    #     [3.0494315263986667 3.464424259522579e-11 41.69589061554819]
    #     [4.054260980052888 0.7624757524487332 53.77764649002646]
    #     [4.054260980052974 -0.7624757523582324 53.77764649002697]
    #     [3.6011693250450887 4.151351781179444e-11 50.039181938658636]
    #     [3.7857414314390576 0.7979344883197637 53.87232980462124]
    #     [3.7857414314395283 -0.7979344882296401 53.872329804622694]
    #     [4.358310593162984 4.520152961124992e-11 53.84035942278267]
    # ]

    # initial init
    k3l.set.mass = 0.9
    k3l.set.l_tether = 50.0
    KiteModels.init_sim!(k3l; prn=false, torque_control=false)
    initial_pos = deepcopy(k3l.pos)
    prn && println("initial_pos")
    for i in 1:3
        @test isapprox(initial_pos[i], [0.0, 0.0, 0.0], atol=tol, rtol=tol)
    end
    for i in 1:3:k3l.num_A
        @test isapprox(initial_pos[i][2], -initial_pos[i+1][2], atol=tol, rtol=tol)
        @test isapprox(initial_pos[i+2][2], 0.0, atol=tol, rtol=tol)
    end
    for i in 4:3:k3l.num_A
        @test initial_pos[i][2] > initial_pos[i-3][2]
    end
    for i in 4:k3l.num_flap_D
        @test initial_pos[i][3] > initial_pos[i-3][3]
    end
    @test initial_pos[k3l.num_A][3] < initial_pos[k3l.num_D][3]
    @test isapprox(norm(initial_pos[k3l.num_E]), k3l.set.l_tether, rtol=0.1)
    @test isapprox(norm(initial_pos[k3l.num_E]), k3l.tether_lengths[3], rtol=0.1)

    if !prn
        # init after changing settings
        k3l.set.mass = 1.0
        k3l.set.l_tether = 51.0
        KiteModels.init_sim!(k3l; prn=false, torque_control=false)
        pos2 = deepcopy(k3l.pos)
        @test isapprox(k3l.tether_lengths[3], 51.0, atol=0.2)
        for i in 4:k3l.num_A
            @test !isapprox(pos2[i], initial_pos[i], atol=tol, rtol=tol)
        end

        # init after changing settings back
        k3l.set.mass = 0.9
        k3l.set.l_tether = 50.0
        KiteModels.init_sim!(k3l; prn=false, torque_control=false)
        pos3 = deepcopy(k3l.pos)
        for i in eachindex(initial_pos)
            @test isapprox(pos3[i], initial_pos[i], atol=tol, rtol=tol)
        end

        # init after changing only initial conditions
        k3l.set.elevation = 84.0
        KiteModels.init_sim!(k3l; prn=false, torque_control=false)
        pos4 = deepcopy(k3l.pos)
        @test isapprox(rad2deg(calc_elevation(k3l)), 84.0, atol=2.0)
        for i in 4:k3l.num_A
            @test !isapprox(pos4[i], initial_pos[i], atol=tol, rtol=tol)
        end

        # init after just stepping
        KiteModels.next_step!(k3l)
        KiteModels.init_sim!(k3l; prn=false, torque_control=false)
        pos5 = deepcopy(k3l.pos)
        for i in eachindex(initial_pos)
            @test isapprox(pos5[i], pos4[i], atol=tol, rtol=tol)
        end
    end

    # TODO: add tests for torque controlled
end

@testset "test_step         " begin
    set_defaults()
    KiteModels.init_sim!(k3l; prn=false, torque_control=false)

    KiteModels.next_step!(k3l)
    # pos2 = [
    #     [2.1851584894939454e-18 -7.505816782042582e-20 4.859625427074995e-18]
    #     [-3.921896020850223e-18 5.946955100867209e-21 6.089875242607995e-19]
    #     [1.0125824674705326e-19 -2.419360367675051e-22 -1.2810021683628663e-18]
    #     [7.788403509281168 0.10572950364523777 6.956685212369045]
    #     [7.788403509183137 -0.1057295036480408 6.956685212478777]
    #     [0.6405221955952705 -1.3689318877296879e-11 8.336480990061558]
    #     [11.611832316831176 0.25053181474596664 16.67400642098399]
    #     [11.61183231692146 -0.2505318147499263 16.674006421019605]
    #     [1.2670834561952096 -2.738412245877073e-11 16.67403178569772]
    #     [12.812598976552383 0.40264305911892356 27.047095510817083]
    #     [12.812598976467324 -0.4026430591286899 27.04709551087292]
    #     [1.8768323853570605 -4.109355859595333e-11 25.01283779368113]
    #     [12.741722921461745 0.5380073055003298 37.4894412805489]
    #     [12.741722921367693 -0.5380073055202815 37.489441280604495]
    #     [2.4679745631752876 -5.4826739072103354e-11 33.35299223068001]
    #     [12.694636474801083 0.6361592296893498 47.93233401589268]
    #     [12.694636474839887 -0.6361592297573209 47.93233401594842]
    #     [3.0390640717513926 -6.858673046731401e-11 41.69455240013717]
    #     [4.039883116055034 0.7624400406478303 53.77613336238382]
    #     [4.039883116056398 -0.7624400408248398 53.77613336238193]
    #     [3.588715003528734 -8.236524620098094e-11 50.03756120467184]
    #     [3.771344262902567 0.7979302697228247 53.870750105023596]
    #     [3.7713442629037255 -0.7979302699002653 53.87075010502097]
    #     [4.343941000937282 -8.824891226050073e-11 53.83912305750433]
    # ]
    # prn && println("pos2")
    for i in 4:k3l.num_A
        @test !isapprox(k3l.pos[i], initial_pos[i], atol = 1e-4)
        # prn ? println(k3l.pos[i]') : @test isapprox(pos2[i,:], k3l.pos[i], atol=tol, rtol=tol)
    end
end

function simulate(steps)
    av_L_C = k3l.get_L_C(k3l.integrator)
    for i in 1:steps
        KiteModels.next_step!(k3l; set_values=[0.0, 0.0, 0.0])
        av_L_C .+= k3l.get_L_C(k3l.integrator)
    end
    av_L_C ./= steps
    return k3l.integrator.iter/steps, av_L_C
end

@testset "test_simulate     " begin
    STEPS = 10
    KiteModels.init_sim!(k3l; prn=false, torque_control=false)
    # println("\nStarting simulation...")
    av_steps, av_L_C = simulate(STEPS)
    prn && println(av_steps)
    if Sys.isapple()
        result = av_steps < 100
        !result && println("isapple, KPS4_3L, steps: $av_steps")
        prn || @test result
    else
        result = av_steps < 100
        !result && println("not apple, KPS4_3L, steps: $av_steps")
        prn || @test result
    end
  
    if prn
        @show k3l.get_L_C(k3l.integrator)
        @show k3l.get_tether_vels(k3l.integrator)
    else
        @test isapprox(av_L_C, [4.205557424653014, 113.23736209049324, 238.29250730249586], rtol=0.1)
        @test isapprox(normalize(k3l.get_L_C(k3l.integrator)) ⋅ normalize(k3l.v_wind), 0.0, atol=0.02)
        @test isapprox(k3l.get_tether_vels(k3l.integrator), [0.0, 0.0, 0.0], atol=0.2)
        @test isapprox(k3l.get_L_C(k3l.integrator)[2], -k3l.get_L_D(k3l.integrator)[2], atol=1e-1)
        @test isapprox(k3l.get_heading(k3l.integrator), 0.0, atol=tol)
    end
    
    # TODO Add testcase with varying reelout speed 
end

# TODO: add testset for sysstate


end

set_data_path(old_path)

nothing