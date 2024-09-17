using REPL.TerminalMenus

options = ["bench = include(\"bench.jl\")",
           "bench_4p = include(\"bench_4p.jl\")",
           "test_init = include(\"test_init.jl\")",
           "plot_cl_cd_plate = include(\"plot_cl_cd_plate.jl\")",
           "compare_kps3_kps4 = include(\"compare_kps3_kps4.jl\")",
           "reel_out_1p = include(\"reel_out_1p.jl\")",
           "reel_out_4p = include(\"reel_out_4p.jl\")",
           "reel_out_4p_torque_control = include(\"reel_out_4p_torque_control.jl\")",
           "simulate_simple = include(\"simulate_simple.jl\")",
           "simulate_steering = include(\"simulate_steering.jl\")",
           "steering_test = include(\"steering_test.jl\")",
           "simple_3l_torque_control = include(\"simple_3l_torque_control.jl\")",
           "quit"]

function menu()
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

menu()