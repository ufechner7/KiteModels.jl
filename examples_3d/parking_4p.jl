using Printf

using Pkg
if ! ("KiteViewers" ∈ keys(Pkg.project().dependencies))
    Pkg.activate("examples_3d")
    pkg"add ControlPlots#main"
    pkg"add KiteUtils"
    pkg"add KiteModels#stable"
end
using KiteModels, KitePodModels, KiteUtils, Rotations, StaticArrays
using ControlPlots, KiteViewers

set = deepcopy(se())

# the following values can be changed to match your interest
dt = 0.05
set.solver="DFBDF"              # IDA or DFBDF
set.linear_solver="GMRES"       # GMRES, LapackDense or Dense
STEPS = 352
PRINT = false
STATISTIC = false
PLOT=true
UPWIND_DIR2       = -pi/2+deg2rad(10)     # Zero is at north; clockwise positive
ZOOM = true
FRONT_VIEW = true
SHOW_KITE = false
# end of user parameter section #

kcu::KCU = KCU(set)
kps4::KPS4 = KPS4(kcu)
viewer::Viewer3D = Viewer3D(SHOW_KITE)

v_time = zeros(STEPS)
v_speed = zeros(STEPS)
v_force = zeros(STEPS)
heading = zeros(STEPS)

# """ 
#     fromKS2EX(vector, orientation)

# Transform a vector (x,y,z) from KiteSensor to Earth Xsens reference frame.

# - orientation in Euler angles (roll, pitch, yaw)
# """
# function fromKS2EX(vector, orientation)
#     roll, pitch, yaw  = orientation[1], orientation[2], orientation[3]
#     rotateYAW = @SMatrix[cos(yaw) -sin(yaw) 0;
#                          sin(yaw)  cos(yaw) 0;
#                              0         0    1]
#     rotatePITCH = @SMatrix[cos(pitch)   0  sin(pitch);
#                              0          1        0;
#                        -sin(pitch)      0  cos(pitch)]
#     rotateROLL = @SMatrix[ 1        0         0;
#                            0   cos(roll) -sin(roll);
#                            0   sin(roll)  cos(roll)]
#     rotateYAW * rotatePITCH * rotateROLL * vector
# end

# """
#     fromEX2EG(vector)

# Transform a vector (x,y,z) from EarthXsens to Earth Groundstation reference frame
# """
# function fromEX2EG(vector)
#     rotateEX2EG = @SMatrix[1  0  0;
#                            0 -1  0;
#                            0  0 -1]
#     rotateEX2EG * vector
# end

# """
#     fromEG2W(vector, down_wind_direction = pi/2.0)

# Transform a vector (x,y,z) from Earth Groundstation to Wind reference frame.
# """
# function fromEG2W2(vector, down_wind_direction = pi/2.0)
#     rotateEG2W =    @SMatrix[cos(down_wind_direction) -sin(down_wind_direction)  0;
#                              sin(down_wind_direction)  cos(down_wind_direction)  0;
#                              0                        0                      1]
#     rotateEG2W * vector
# end

# function calc_heading_w2(orientation, down_wind_direction = pi/2.0)
#     # create a unit heading vector in the xsense reference frame
#     heading_sensor =  SVector(1, 0, 0)
#     # rotate headingSensor to the Earth Xsens reference frame
#     headingEX = fromKS2EX(heading_sensor, orientation)
#     # rotate headingEX to earth groundstation reference frame
#     headingEG = fromEX2EG(headingEX)
#     # rotate headingEG to headingW and convert to 2d HeadingW vector
#     fromEG2W2(headingEG, down_wind_direction)
# end

# """
#     calc_heading(orientation, elevation, azimuth; upwind_dir=-pi/2, respos=true)

# Calculate the heading angle of the kite in radians. The heading is the direction
# the nose of the kite is pointing to. 
# If respos is true the heading angle is defined in the range of 0 .. 2π,
# otherwise in the range -π .. π
# """
# function calc_heading2(orientation, elevation, azimuth; upwind_dir=-pi/2, respos=true)
#     down_wind_direction = wrap2pi(upwind_dir + π)
#     headingSE = fromW2SE(calc_heading_w2(orientation, down_wind_direction), elevation, azimuth)
#     angle = atan(headingSE.y, headingSE.x) # - π
#     if angle < 0 && respos
#         angle += 2π
#     end
#     angle
# end

# function calc_heading2(s::KPS4; upwind_dir=upwind_dir(s))
#     orientation = orient_euler(s)
#     elevation = calc_elevation(s)
#     azimuth = calc_azimuth(s)
#     println("azimuth: ", rad2deg(azimuth))
#     calc_heading2(orientation, elevation, azimuth; upwind_dir)
# end

function simulate(integrator, steps, plot=true)
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
            steering = 0.05
        end
        set_depower_steering(kps4.kcu, kps4.depower, steering)

        KiteModels.next_step!(kps4, integrator; set_speed, dt, upwind_dir=UPWIND_DIR2)
        iter += kps4.iter
        if plot
            reltime = i*dt-dt
            if mod(i, 5) == 1
                plot2d(kps4.pos, reltime; zoom=true, front=FRONT_VIEW, 
                       segments=set.segments, fig="front_view") 
                sleep(0.05)           
            end
        end
        sys_state = SysState(kps4)
        q = QuatRotation(sys_state.orient)
        q_viewer = AngleAxis(-π/2, 0, 1, 0) * q
        sys_state.orient .= Rotations.params(q_viewer)
        KiteViewers.update_system(viewer, sys_state; scale = 0.08, kite_scale=3)
        # if wrap2pi(calc_heading(kps4)) > 0 && i > 100; break; end
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
lift, drag = KiteModels.lift_drag(kps4)
println("lift, drag  [N]: $(round(lift, digits=2)), $(round(drag, digits=2))")

println("v_wind: $(kps4.v_wind)")
println("UPWIND_DIR2: $(rad2deg(UPWIND_DIR2))°")
pos = pos_kite(kps4)
println("pos_y: $(round(pos[2], digits=2))")
# for an  UPWIND_DIR2 of -80°, pos_y must be negative, also v_wind[2] must be negative
# this is OK

# print heading
println("heading: $(round(heading[STEPS], digits=2))°")
plot(v_time, heading; xlabel="time [s]", ylabel="heading [°]", fig="heading")