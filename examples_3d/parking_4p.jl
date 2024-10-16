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

function euler2rot(roll, pitch, yaw)
    φ      = roll
    R_x = [1    0       0;
              0  cos(φ) -sin(φ);
              0  sin(φ)  cos(φ)]
    θ      = pitch          
    R_y = [ cos(θ)  0  sin(θ);
                 0     1     0;
              -sin(θ)  0  cos(θ)]
    ψ      = yaw
    R_z = [cos(ψ) -sin(ψ) 0;
              sin(ψ)  cos(ψ) 0;
                 0       0   1]
    R   = R_z * R_y * R_x
    return R
end

"""
    quat2viewer(q::QuatRotation)
    quat2viewer(rot::AbstractMatrix)
    quat2viewer(orient::AbstractVector)

Convert the quaternion q to the viewer reference frame. It can also be passed
as a rotation matrix or as 4-element vector [w,i,j,k], where w is the real part
and i, j, k are the imaginary parts of the quaternion.
"""
quat2viewer_(rot::AbstractMatrix) = quat2viewer(QuatRotation(rot))
quat2viewer_(orient::AbstractVector) = quat2viewer(QuatRotation(orient))
function quat2viewer_(q::QuatRotation)
    # 1. get reference frame
    rot = inv(RotMatrix{3}(q))
    x = enu2ned(rot[1,:])
    y = enu2ned(rot[2,:])
    z = enu2ned(rot[3,:])
    # 2. convert it using the old method
    ax = [0, 1, 0] # in ENU reference frame this is pointing to the south
    ay = [1, 0, 0] # in ENU reference frame this is pointing to the west
    az = [0, 0, -1] # in ENU reference frame this is pointing down
    rotation = rot3d(ax, ay, az, x, y, z)
    q_old = QuatRotation(rotation)
    x = [0,  1.0, 0]
    y = [1.0,  0, 0]
    z = [0,    0, -1.0]
    x, y, z = q_old*x, q_old*y, q_old*z
    rot = calc_orient_rot(x, y, z; viewer=true, ENU=false)
    q = QuatRotation(rot)
    return Rotations.params(q)
end

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

        sys_state.orient = quat2viewer_(q)
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

# # print heading
# println("heading: $(round(heading[STEPS], digits=2))°")
# plot(v_time, heading; xlabel="time [s]", ylabel="heading [°]", fig="heading")