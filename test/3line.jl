using KiteModels
using ControlPlots
using LinearAlgebra, StaticArrays
using REPL

update_settings()
kcu = KCU(se())
s = KPS4_3L(kcu)
integrator = KiteModels.init_sim!(s, stiffness_factor=0.04, prn=true, integrator_history=nothing)

# Run listener as separate task using channels, put keypresses in channel for main loop
function key_listener(c::Channel)
    t = REPL.TerminalMenus.terminal
    while true
        REPL.Terminals.raw!(t, true) || error("unable to switch to raw mode")
        keypress = Char(REPL.TerminalMenus.readkey(t.in_stream))
        REPL.Terminals.raw!(t, false) || error("unable to switch back from raw mode")
        put!(c, keypress)
    end
end

function main()
    channel = Channel(key_listener, 10) # Start task, 10 is buffer size for channel
    stop = false
    v_ro=[0.0, 0.0, 0.0]
    t = 0
    dt = 0.05
    flag = true
    last_time = 0.0
    while !stop
        if flag
            t += dt
            @time KiteModels.next_step!(s, integrator, v_ro=v_ro, dt=dt)
            plot2d(s.pos, t; zoom=false, front=false, segments=s.set.segments)
            println(v_ro)
            println("δ_left ", s.δ_left)
            println("δ_right ", s.δ_right)
            println("tether lengths ", s.tether_lengths)
            println("heading ", calc_heading(s))
            while time() - last_time < dt
                sleep(0.01)
            end
            last_time = time()
        end
    
        while !isempty(channel)
            c = take!(channel)
            v_ro = [0.0, 0.0, 0.0]
            if c == 'z'
                stop = true
                close(channel)
                break
            else
                flag=true
                if c == 'a'
                    v_ro[2] = 0.5
                end
                if c == 's'
                    v_ro[2] = 0.5
                    v_ro[3] = 0.5
                end
                if c == 'd'
                    v_ro[3] = 0.5
                end
                if c == 'q'
                    println("reel out left")
                    v_ro[2] = -0.5
                end
                if c == 'w'
                    v_ro[2] = -0.5
                    v_ro[3] = -0.5
                end
                if c == 'e'
                    println("reel out right")
                    v_ro[3] = -0.5
                end
                break
            end
        end
    end
end

main()

# i=0
# while true # steer kite
#     global i += 1
#     dt = 0.2
#     plot2d(s.pos, i*dt; zoom=false, front=true, segments=s.set.segments)
#     println(i)
#     println("δ_left ", s.δ_left)
#     println("δ_right ", s.δ_right)
#     println("lengths ", s.tether_lengths)
#     # println("Pos \t", s.pos)
#     println(check_key_pressed())
#     if i < 70
#         @time KiteModels.next_step!(s, integrator, v_ro=[0.0,0.0,-0.5], dt=dt)
#     elseif i < 140
#         @time KiteModels.next_step!(s, integrator, v_ro=[0.0,-0.5,0.0], dt=dt)
#     else
#         @time KiteModels.next_step!(s, integrator, v_ro=[0.0,0.0,0.0], dt=dt)
#     end
# end