# Copyright (c) 2022, 2024 Uwe Fechner
# SPDX-License-Identifier: MIT
using REPL.TerminalMenus

options = ["bench = include(\"bench.jl\")",
           "bench_4p = include(\"bench_4p.jl\")",
           "test_init_1p = include(\"test_init_1p.jl\")",
           "test_init_4p = include(\"test_init_4p.jl\")",
           "plot_cl_cd_plate = include(\"plot_cl_cd_plate.jl\")",
           "plot_side_cl = include(\"plot_side_cl.jl\")",
           "compare_kps3_kps4 = include(\"compare_kps3_kps4.jl\")",
           "reel_out_1p = include(\"reel_out_1p.jl\")",
           "reel_out_4p = include(\"reel_out_4p.jl\")",
           "reel_out_4p_torque_control = include(\"reel_out_4p_torque_control.jl\")",
           "simulate_simple = include(\"simulate_simple.jl\")",
           "simulate_steering = include(\"simulate_steering.jl\")",
           "steering_test_1p = include(\"steering_test_1p.jl\")",
           "steering_test_4p = include(\"steering_test_4p.jl\")",
           "ram_air_kite = SIMPLE=false; include(\"ram_air_kite.jl\")",
           "simple_ram_air_kite = SIMPLE=true; include(\"ram_air_kite.jl\")",
           "lin_ram_model = include(\"lin_ram_model.jl\")",
           "calc_spectrum = include(\"calc_spectrum.jl\")",
           "plot_spectrum_ = include(\"plot_spectrum.jl\")",
           "calculate_rotational_inertia = include(\"calculate_rotational_inertia.jl\")",
           "quit"]

function example_menu()
    active = true
    while active
        menu = RadioMenu(options, pagesize=8)
        choice = request("\nChoose function to execute or `q` to quit: ", menu)

        if choice != -1 && choice != length(options)
            eval(Meta.parse(options[choice]))
        else
            println("Left menu. Press <ctrl><d> to quit Julia!")
            active = false
        end
    end
end

example_menu()