# TODO: calculate the rotation between q and q2
# - first add the function to calculate the rotation between two quaternions

using Printf

using Pkg, Timers
tic()
if ! ("KiteViewers" ∈ keys(Pkg.project().dependencies))
    Pkg.activate("examples_3d")
    pkg"add KiteUtils#main"
    pkg"add KiteModels#main"
end
using KiteModels, KitePodModels, KiteUtils, Rotations, StaticArrays
using ControlPlots, KiteViewers
toc()

set = deepcopy(se())

# the following values can be changed to match your interest
dt = 0.05
set.solver="DFBDF"              # IDA or DFBDF
set.linear_solver="GMRES"       # GMRES, LapackDense or Dense
STEPS = 352
PRINT = false
STATISTIC = false
PLOT=false
UPWIND_DIR2       = -pi/2+deg2rad(10)     # Zero is at north; clockwise positive
ZOOM = true
FRONT_VIEW = true
SHOW_KITE = true
# end of user parameter section #

kcu::KCU = KCU(set)
kps4::KPS4 = KPS4(kcu)
viewer::Viewer3D = Viewer3D(SHOW_KITE)

v_time = zeros(STEPS)
v_speed = zeros(STEPS)
v_force = zeros(STEPS)
heading = zeros(STEPS)

function simulate(integrator, steps, plot=PLOT)
    iter = 0
    for i in 1:steps
        acc = 0.0
        v_time[i] = kps4.t_0
        v_speed[i] = kps4.v_reel_out
        v_force[i] = winch_force(kps4)
        heading[i] = rad2deg(wrap2pi(calc_heading(kps4)))
        set_speed = kps4.sync_speed+acc*dt
        if PRINT
            lift, drag = KiteModels.lift_drag(kps4)
            @printf "%.2f: " round(integrator.t, digits=2)
            println("lift, drag  [N]: $(round(lift, digits=2)), $(round(drag, digits=2))")
        end
        steering = 0
        if i > 200
            steering = -0.0482
        end
        set_depower_steering(kps4.kcu, kps4.depower, steering)

        KiteModels.next_step!(kps4, integrator; set_speed, dt, upwind_dir=UPWIND_DIR2)
        iter += kps4.iter
        reltime = i*dt-dt
        if mod(i, 5) == 1
            if plot
                plot2d(kps4.pos, reltime; zoom=true, front=FRONT_VIEW, 
                    segments=set.segments, fig="front_view") 
            end
            sleep(0.05)           
        end
        sys_state = SysState(kps4)
        q = QuatRotation(sys_state.orient)
        # roll, pitch, yaw = quat2euler(q)
        # println("Yaw: ", rad2deg(yaw), ", Pitch: ", rad2deg(pitch), ", Roll: ", rad2deg(roll))

        sys_state.orient = quat2viewer(q)
        KiteViewers.update_system(viewer, sys_state; scale = 0.08, kite_scale=3)
    end
    iter / steps
end

integrator = KiteModels.init_sim!(kps4, delta=0, stiffness_factor=0.5, prn=STATISTIC)

println("\nStarting simulation...")
simulate(integrator, STEPS)
if PLOT
    p = plotx(v_time[1:STEPS-100], v_speed[1:STEPS-100], v_force[1:STEPS-100]; ylabels=["v_reelout  [m/s]","tether_force [N]"], fig="winch")
    # display(p)
end
# lift, drag = KiteModels.lift_drag(kps4)
# println("lift, drag  [N]: $(round(lift, digits=2)), $(round(drag, digits=2))")

# println("v_wind: $(kps4.v_wind)")
# println("UPWIND_DIR2: $(rad2deg(UPWIND_DIR2))°")
# pos = pos_kite(kps4)
# println("pos_y: $(round(pos[2], digits=2))")
# # for an  UPWIND_DIR2 of -80°, pos_y must be negative, also v_wind[2] must be negative
# # this is OK

# print heading
println("heading: $(round(heading[STEPS], digits=2))°")
plot(v_time, heading; xlabel="time [s]", ylabel="heading [°]", fig="heading")