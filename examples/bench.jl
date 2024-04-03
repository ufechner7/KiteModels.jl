using Printf
using KiteModels, KitePodModels, KiteUtils

if ! @isdefined kcu;  const kcu = KCU(se());   end
if ! @isdefined kps4; const kps4 = KPS4(kcu); end

# the following values can be changed to match your interest
dt = 0.05
STEPS = 600
PLOT = false
FRONT_VIEW = false
ZOOM = true
PRINT = false
STATISTIC = false
# end of user parameter section #

function simulate(integrator, steps, plot=false)
    start = integrator.p.iter
    for i in 1:steps
        if PRINT
            lift, drag = KiteModels.lift_drag(kps4)
            @printf "%.2f: " round(integrator.t, digits=2)
            println("lift, drag  [N]: $(round(lift, digits=2)), $(round(drag, digits=2))")
        end

        KiteModels.next_step!(kps4, integrator, dt=dt)
    end
    (integrator.p.iter - start) / steps
end

integrator = KiteModels.init_sim!(kps4, stiffness_factor=0.04, prn=STATISTIC)

println("\nStarting simulation...")
simulate(integrator, 100)
runtime = @elapsed av_steps = simulate(integrator, STEPS-100)
println("\nTotal simulation time: $(round(runtime, digits=3)) s")
speed = (STEPS-100) / runtime * dt
println("Simulation speed: $(round(speed, digits=2)) times realtime.")
lift, drag = KiteModels.lift_drag(kps4)
println("lift, drag  [N]: $(round(lift, digits=2)), $(round(drag, digits=2))")
println("Average number of callbacks per time step: $av_steps")

# Ryzen 7950X, GMRES solver
# Total simulation time: 0.712 s
# Simulation speed: 35.11 times realtime.
# lift, drag  [N]: 597.4, 129.32
# Average number of callbacks per time step: 378.918

# Ryzen 7950X, LapackDense solver
# Total simulation time: 2.216 s
# Simulation speed: 11.28 times realtime.
# lift, drag  [N]: 597.56, 129.31
# Average number of callbacks per time step: 1413.298

# Total simulation time: 0.332 s
# Simulation speed: 75.25 times realtime.
# lift, drag  [N]: 755.2, 168.82
# Average number of callbacks per time step: 227.5

# Total simulation time: 0.101 s
# Simulation speed: 248.17 times realtime.
# lift, drag  [N]: 811.03, 183.15
# Average number of callbacks per time step: 68.238