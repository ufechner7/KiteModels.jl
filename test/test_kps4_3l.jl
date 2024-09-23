using Test, BenchmarkTools, StaticArrays, LinearAlgebra, KiteUtils
using KiteModels, KitePodModels

kcu_3l::KCU = KCU(se("system_3l.yaml"))
kcu_3l.set.winch_model = "AsyncMachine"
kps4_3l::KPS4_3L = KPS4_3L(kcu_3l)

pos, vel = nothing, nothing

@testset verbose = true "KPS4_3L tests...." begin

tol::Float64 = 1e-7

function set_defaults()
    kps4_3l.set = update_settings()
    # kps4_3l.set.abs_tol = tol
    # kps4_3l.set.rel_tol = tol
    KiteModels.clear!(kps4_3l)
end

set_defaults()

@testset "test_init         " begin
    set_defaults()
    [kps4_3l.pos[i] .= 0 for i in 1:kps4_3l.num_A]

    initial_pos = [
        [-6.696548532834703e-17 4.3241435580763785e-18 9.697253528187506e-17]
        [5.400221440010069e-17 -5.599654834593977e-18 -1.6989080336075647e-16]
        [1.7142142314719172e-16 1.1031542481779438e-17 -3.843686938071738e-17]
        [7.734084348247082 0.10645077963181254 7.0170346494179]
        [7.73408434810422 -0.10645077960863832 7.017034649575712]
        [0.5529745919321083 1.0758628300353887e-12 8.34402114746499]
        [11.529524989680493 0.2512157975094456 16.745335844195292]
        [11.52952498939103 -0.25121579750242073 16.745335844410047]
        [1.0925239738447141 2.1507195189834437e-12 16.688930158548825]
        [12.710900384978052 0.4028043436965661 27.120672081654764]
        [12.710900385218649 -0.40280434369885054 27.12067208180902]
        [1.6159339099375558 3.2228808426197586e-12 25.034876042232654]
        [12.55800143136617 0.5359264415265507 37.562184313595196]
        [12.558001431333789 -0.5359264415322457 37.5621843137454]
        [2.1215309769617363 4.288088147188535e-12 33.38192903034019]
        [12.163348198024075 0.6315024573573741 47.997768282326504]
        [12.163348198128432 -0.6315024573607547 47.99776828248185]
        [2.6080126625657205 5.33665863009136e-12 41.730127011347996]
        [3.485808990663181 0.7561165912989308 53.807750977747794]
        [3.4858089906635255 -0.7561165912855672 53.807750977748114]
        [3.074276396486973 6.348525440467502e-12 50.07948775187569]
        [3.2227217569476525 0.7979433702173507 53.914337224614926]
        [3.2227217569479163 -0.7979433702039685 53.91433722461505]
        [3.7955460652648867 6.784280233584585e-12 53.88763482966672]
    ]

    # initial init
    kps4_3l.set.mass = 0.9
    kps4_3l.set.l_tether = 50.0
    KiteModels.init_sim!(kps4_3l; prn=true, torque_control=false)
    pos1 = deepcopy(kps4_3l.pos)
    for i in eachindex(pos1)
        # println(pos1[i]')
        @test isapprox(pos1[i], initial_pos[i, :], atol=tol, rtol=tol)
    end

    # init after changing settings
    kps4_3l.set.mass = 1.0
    kps4_3l.set.l_tether = 51.0
    KiteModels.init_sim!(kps4_3l; prn=true, torque_control=false)
    pos2 = deepcopy(kps4_3l.pos)
    @test isapprox(kps4_3l.tether_lengths[3], 51.0, atol=0.1)
    for i in 4:kps4_3l.num_A
        @test !isapprox(pos2[i], initial_pos[i, :], atol=tol, rtol=tol)
    end

    # init after changing settings back
    kps4_3l.set.mass = 0.9
    kps4_3l.set.l_tether = 50.0
    KiteModels.init_sim!(kps4_3l; prn=true, torque_control=false)
    pos3 = deepcopy(kps4_3l.pos)
    for i in eachindex(pos1)
        @test isapprox(pos3[i], initial_pos[i, :], atol=tol, rtol=tol)
    end

    # init after changing only initial conditions
    kps4_3l.set.elevation = 84.0
    KiteModels.init_sim!(kps4_3l; prn=true, torque_control=false)
    pos4 = deepcopy(kps4_3l.pos)
    @test isapprox(rad2deg(calc_elevation(kps4_3l)), 84.0, atol=2.0)
    for i in 4:kps4_3l.num_A
        @test !isapprox(pos4[i], initial_pos[i, :], atol=tol, rtol=tol)
    end

    # init after just stepping
    KiteModels.next_step!(kps4_3l)
    KiteModels.init_sim!(kps4_3l; prn=true, torque_control=false)
    pos5 = deepcopy(kps4_3l.pos)
    for i in eachindex(pos1)
        @test isapprox(pos5[i], pos4[i], atol=tol, rtol=tol)
    end

    # TODO: add tests for torque controlled
end

@testset "test_step         " begin
    set_defaults()
    KiteModels.init_sim!(kps4_3l; prn=true, torque_control=false)

    KiteModels.next_step!(kps4_3l)
    pos2 = [
        [7.010919237850931e-17 -1.758540841573567e-18 -3.106375645731575e-16]
        [2.416999885990368e-16 -5.136990049706791e-18 -4.573687906382295e-16]
        [2.1260442326447698e-16 3.79227781182154e-18 -1.2605807013224936e-17]
        [7.7366648156385285 0.10648445927453026 7.014172016683147]
        [7.736664815487313 -0.10648445925147394 7.014172016849973]
        [0.5515251858606224 4.5042497465162454e-12 8.343730721771747]
        [11.549988820790695 0.2511510842059863 16.73546858188156]
        [11.549988820495258 -0.2511510841985439 16.73546858210443]
        [1.0895990335326886 9.003206820429557e-12 16.688348647382867]
        [12.74363989706399 0.40266958155569554 27.109391705347388]
        [12.743639897298973 -0.4026695815564375 27.109391705508912]
        [1.611547863591678 1.3475755138014307e-11 25.03399964161336]
        [12.587533245051867 0.5358335326506541 37.55084457629541]
        [12.587533245022149 -0.5358335326518386 37.5508445764528]
        [2.115732734116989 1.7854934990431922e-11 33.38075158452793]
        [12.145268940696658 0.6318368320078857 47.98450045887626]
        [12.145268940806899 -0.6318368319977479 47.98450045903954]
        [2.600864955214557 2.2119111459596548e-11 41.72864152420909]
        [3.475528597445595 0.7561065325592121 53.806087359523104]
        [3.4755285974458365 -0.7561065325014515 53.80608735952352]
        [3.0657203433763534 2.7647139238117296e-11 50.07769408442735]
        [3.2123980606508624 0.7979394009595864 53.912564266538254]
        [3.2123980606510965 -0.7979394009018203 53.912564266538645]
        [3.7852451249488954 2.8961939866696227e-11 53.88617298127155]
    ]
    for i in eachindex(kps4_3l.pos)
        println(kps4_3l.pos[i]')
        # @test isapprox(pos2[i,:], kps4_3l.pos[i], atol=tol, rtol=tol)
    end
    println(kps4_3l.L_C)
    # @test isapprox(kps4_3l.L_C, [-0.9050171285048465, 146.0015097898251, 307.6023126186097], atol=tol, rtol=tol)
    @test isapprox(normalize(kps4_3l.L_C) ⋅ normalize(kps4_3l.v_wind), 0.0, atol=1e-2)
end

function simulate(steps)
    for i in 1:steps
        KiteModels.next_step!(kps4_3l; set_values=[0.0, 0.0, 0.0])
    end
    return kps4_3l.integrator.iter/steps
end

@testset "test_simulate     " begin
    STEPS = 20
    KiteModels.init_sim!(kps4_3l; prn=true, torque_control=false)
    # println("\nStarting simulation...")
    av_steps = simulate(STEPS)
    # println(av_steps)
    if Sys.isapple()
        println("isapple $av_steps")
        @test isapprox(av_steps, 3.65, atol=1.0)
    else
        println("not apple $av_steps")
        @test isapprox(av_steps, 3.65, atol=1.0)
    end
  
    # @show kps4_3l.L_C
    # @show kps4_3l.reel_out_speeds
    @test isapprox(kps4_3l.L_C, [0.5481297824282668, 147.23411730075136, 310.2882790830033], atol=1.0)
    @test isapprox(kps4_3l.reel_out_speeds, [0.0, 0.0, 0.0], atol=tol)
    @test isapprox(kps4_3l.L_C[2], -kps4_3l.L_D[2], atol=1e-3)
    
    # TODO Add testcase with varying reelout speed 
end

# TODO: add testset for sysstate

end
nothing