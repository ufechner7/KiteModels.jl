using REPL.TerminalMenus

options = ["plot_cl_cd_plate = include(\"plot_cl_cd_plate.jl\")",
           "plot_side_cl = include(\"plot_side_cl.jl\")",
           "steering_test_4p = include(\"steering_test_4p.jl\")",
           "plot_parking_test = include(\"plot_parking_test.jl\")",
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