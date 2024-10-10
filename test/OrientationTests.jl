using Pkg

using KiteModels
using KitePodModels
using KiteUtils
using LinearAlgebra
using StaticArrays
using Test
using Printf

function create_kite_model(x, y, z, pos)
    kcu::KCU = KCU(se())
    s::KPS4 = KPS4(kcu)

    s.x = x
    s.y = y
    s.z = z

    s.pos[end-2][begin] = pos[begin]
    s.pos[end-2][begin+1] = pos[begin+1]
    s.pos[end-2][begin+2] = pos[begin+2]
    set_v_wind_ground!(s)
    s
end  

function create_kite_model(x, y, z, pos, wind)
    kcu::KCU = KCU(se())
    s::KPS4 = KPS4(kcu)

    s.x = x
    s.y = y
    s.z = z

    set_v_wind_ground!(s, pos[begin+2], wind_dir=deg2rad(wind))

    s.pos[end-2][begin] = pos[begin]
    s.pos[end-2][begin+1] = pos[begin+1]
    s.pos[end-2][begin+2] = pos[begin+2]
    
    s
end  

function obtain_results(s)
    roll, pitch, yaw = orient_euler(s)
    elevation = calc_elevation(s)
    azimuth = calc_azimuth(s)
    heading = calc_heading(s)
    return roll, pitch, yaw, azimuth, elevation, heading
end

function print_results(test, roll, pitch, yaw, elevation, azimuth, heading)
    print(test, "\n")
    @printf("Roll: %0.2f, Pitch: %0.2f, Yaw: %0.2f\n", rad2deg(roll), rad2deg(pitch), rad2deg(yaw))
    @printf("Elevation: %0.2f, Azimuth: %0.2f, Heading: %0.2f\n", rad2deg(elevation), rad2deg(azimuth), rad2deg(heading))
    
end

function test_attitude(roll, pitch, yaw, correct)
    c_roll, c_pitch, c_yaw = correct
    return rad2deg(roll) == c_roll && rad2deg(pitch) == c_pitch && rad2deg(yaw) == c_yaw
end
    
function test_SE_attitude(azimuth, elevation, heading, correct)
    c_azimuth, c_elevation, c_heading = correct
    return rad2deg(azimuth) = c_azimuth && rad2deg(elevation) == c_elevation && rad2deg(heading) == c_heading
end

"The Heading here is hard to determine, so we should either egree that it defaults to a value of 0 or that it throws an execption`"
# function test_zero_euler()
#     s = create_kite_model((0, 1, 0), (1, 0, 0), (0, 0, -1), (0, 1,0))
#     roll, pitch, yaw, azimuth, elevation, heading = obtain_results(s)
#     print_results(test_zero_euler, roll, pitch, yaw, elevation, azimuth, heading)
#     return roll == pitch  == yaw == elevation == azimuth == 0 && heading == 0 
           
# end 

function test_elevation_45()
    s = create_kite_model((0, 1, 0), (1, 0, 0), (0, 0, -1), # Orientation
                          (0, sqrt(2)/2, sqrt(2)/2))        # Pos
    roll, pitch, yaw, azimuth, elevation, heading = obtain_results(s)
    print_results(test_elevation_45, roll, pitch, yaw, elevation, azimuth, heading)
    return test_attitude(roll, pitch, yaw, (0, 0, 0)) && test_SE_attitude(0, 45, 180)
end
    
function test_elevation_60()
    s = create_kite_model((0, 1, 0), (1, 0, 0), (0, 0, -1), # Orientation
                          (0, 0.5, sqrt(3)/2))              # Pos
    roll, pitch, yaw, azimuth, elevation, heading = obtain_results(s)
    print_results(test_elevation_60, roll, pitch, yaw, elevation, azimuth, heading)
    return test_attitude(roll, pitch, yaw, (0, 0, 0)) && test_SE_attitude(0, 60, 180)
end

function test_elevation_minus_60()
    s = create_kite_model((0, 1, 0), (1, 0, 0), (0, 0, -1), # Orientation
                          (0, 0.5, -sqrt(3)/2))             # Pos
    roll, pitch, yaw, azimuth, elevation, heading = obtain_results(s)
    print_results(test_elevation_minus_60, roll, pitch, yaw, elevation, azimuth, heading)
    return test_attitude(roll, pitch, yaw, (0, 0, 0)) && test_SE_attitude(0, -60, 0)
end

function test_azimuth_45()
    s = create_kite_model((0, 1, 0), (1, 0, 0), (0, 0, -1), # Orientation
                          (-sqrt(2)/2, sqrt(2)/2, 0))       # Pos
    roll, pitch, yaw, azimuth, elevation, heading = obtain_results(s)
    print_results(test_azimuth_45, roll, pitch, yaw, elevation, azimuth, heading)
    return test_attitude(roll, pitch, yaw, (0, 0, 0)) && test_SE_attitude(45, 0, 90)
end

function test_azimuth_60()
    s = create_kite_model((0, 1, 0), (1, 0, 0), (0, 0, -1), # Orientation
                          (-sqrt(3)/2, 1/2, 0))             # Pos
    roll, pitch, yaw, azimuth, elevation, heading = obtain_results(s)
    print_results(test_azimuth_60, roll, pitch, yaw, elevation, azimuth, heading)
    return test_attitude(roll, pitch, yaw, (0, 0, 0)) && test_SE_attitude(60, 0, 90)
end

function test_azimuth_minus_60()
    s = create_kite_model((0, 1, 0), (1, 0, 0), (0, 0, -1), # Orientation
                           (sqrt(3)/2, 1/2, 0))             # Pos
    roll, pitch, yaw, azimuth, elevation, heading = obtain_results(s)
    print_results(test_azimuth_60, roll, pitch, yaw, elevation, azimuth, heading)
    return test_attitude(roll, pitch, yaw, (0, 0, 0)) && test_SE_attitude(-60, 0, 270)
end

function test_pitch_azimuth()
    s = create_kite_model((0, sqrt(2)/2, -sqrt(2)/2), (1, 0, 0), (0, -sqrt(2)/2, -sqrt(2)/2),   # Orientation
                          (-sqrt(3)/2, 1/2, 0))                                                 # Pos
    roll, pitch, yaw, azimuth, elevation, heading = obtain_results(s)
    print_results(test_pitch_azimuth, roll, pitch, yaw, elevation, azimuth, heading)
    return test_attitude(roll, pitch, yaw, (0, 45, 0)) && test_SE_attitude(60, 0, 315)
end

function test_yaw_elevation()
    s = create_kite_model((-sqrt(2)/2, sqrt(2)/2, 0), (sqrt(2)/2, sqrt(2)/2, 0), (0, 0, -1),    # Orientation
                          (0, 1/2, sqrt(3)/2))                                                  # Pos
    roll, pitch, yaw, azimuth, elevation, heading = obtain_results(s)
    print_results(test_yaw_elevation, roll, pitch, yaw, elevation, azimuth, heading)
    return test_attitude(roll, pitch, yaw, (0, 0, 45)) && test_SE_attitude(0, 60, 135)
end

function test_azimuth_45_elevation_45()
    s = create_kite_model((0, 1, 0), (1, 0, 0), (0, 0, -1), # Orientation
                          (-1/2, 1/2, sqrt(2)/2))           # Pos
    roll, pitch, yaw, azimuth, elevation, heading = obtain_results(s)
    print_results(test_yaw_elevation, roll, pitch, yaw, elevation, azimuth, heading)
    return test_attitude(roll, pitch, yaw, (0, 0, 0)) && test_SE_attitude(0, 60, 135)
end

function test_roll_azimuth_elevation()
    s = create_kite_model((0, 1, 0), (0, 0, 1), (1, 0, 0), # Orientation
                          (-1/2, 1/2, sqrt(2)/2))           # Pos
    roll, pitch, yaw, azimuth, elevation, heading = obtain_results(s)
    print_results(test_yaw_elevation, roll, pitch, yaw, elevation, azimuth, heading)
    return test_attitude(roll, pitch, yaw, (0, 90, 0)) && test_SE_attitude(0, 60, 135)
end

function test_base_orientation()
    s = create_kite_model((0, 0, 1), (-1, 0, 0), (0, -1, 0), # Orientation
                          (0, 1, 0))                               # Pos
    roll, pitch, yaw, azimuth, elevation, heading = obtain_results(s)
    print_results(test_yaw_elevation, roll, pitch, yaw, elevation, azimuth, heading)
    return test_attitude(roll, pitch, yaw, (0, -90, 180)) && test_SE_attitude(0, 0, 0)
end

function test_base_orientation_heading_90()
    s = create_kite_model((-1, 0, 1), (0, 0, -1), (0, -1, 0), # Orientation
                          (0, 1, 0))                          # Pos
    roll, pitch, yaw, azimuth, elevation, heading = obtain_results(s)
    print_results(test_yaw_elevation, roll, pitch, yaw, elevation, azimuth, heading)
    return test_attitude(roll, pitch, yaw, (-90, 0, 180)) && test_SE_attitude(0, 0, 90)
end

function test_wind_dir_45()
    s = create_kite_model((-1, 0, 1), (0, 0, -1), (0, -1, 0), # Orientation
                          (0, 1, 0),                          # Pos
                          45)                                 # Wind dir  
    roll, pitch, yaw, azimuth, elevation, heading = obtain_results(s)
    print_results(test_yaw_elevation, roll, pitch, yaw, elevation, azimuth, heading)
    return test_attitude(roll, pitch, yaw, (0, 0, 0)) && test_SE_attitude(0, -45, 0)
end

function test_wind_dir_60()
    s = create_kite_model((-1, 0, 1), (0, 0, -1), (0, -1, 0), # Orientation
                          (0, 1, 0),                          # Pos
                          60)                                 # Wind dir  
    roll, pitch, yaw, azimuth, elevation, heading = obtain_results(s)
    print_results(test_yaw_elevation, roll, pitch, yaw, elevation, azimuth, heading)
    return test_attitude(roll, pitch, yaw, (0, 0, 0)) && test_SE_attitude(0, -60, 0)
end

function test_wind_dir_minus_60()
    s = create_kite_model((-1, 0, 1), (0, 0, -1), (0, -1, 0), # Orientation
                          (0, 1, 0),                          # Pos
                          -60)                                 # Wind dir  
    roll, pitch, yaw, azimuth, elevation, heading = obtain_results(s)
    print_results(test_yaw_elevation, roll, pitch, yaw, elevation, azimuth, heading)
    return test_attitude(roll, pitch, yaw, (0, 0, 0)) && test_SE_attitude(0, 60, 0)
end

@test test_elevation_45()
@test test_elevation_60()
@test test_elevation_minus_60()
@test test_azimuth_45()
@test test_azimuth_60()
@test test_azimuth_minus_60()
@test test_pitch_azimuth()
@test test_yaw_elevation()
@test test_azimuth_45_elevation_45()
@test test_roll_azimuth_elevation()
@test test_base_orientation()
@test test_base_orientation_heading_90()
@test test_wind_dir_45()
@test test_wind_dir_60()
@test test_wind_dir_minus_60()