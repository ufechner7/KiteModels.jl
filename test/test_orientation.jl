# SPDX-FileCopyrightText: 2025 Uwe Fechner, Daan van Wolffelaar
# SPDX-License-Identifier: MIT

using KiteUtils
using KiteModels
using KitePodModels

using LinearAlgebra
using StaticArrays
using Test

"""
    create_kite_model(x, y, z, pos)

Create a kite model with a given kite reference frame and kite position.

x, y, z:    Kite reference frame in ENU coordinates
pos:        Kite position in ENU coordinates
"""
function create_kite_model(x, y, z, pos)
    kcu::KCU = KCU(se())
    s::KPS4 = KPS4(kcu)

    s.x = x
    s.y = y
    s.z = z

    KiteModels.set_v_wind_ground!(s, pos[begin+2], deg2rad(-90))

    s.pos[end-2][begin] = pos[begin]
    s.pos[end-2][begin+1] = pos[begin+1]
    s.pos[end-2][begin+2] = pos[begin+2]
    s
end  

"""
    create_kite_model(x, y, z, pos, upwind_dir_deg)

Create a kite model with a given kite reference frame, kite position and
upwind direction.

x, y, z:        Kite reference frame in ENU coordinates
pos:            Kite position in ENU coordinates
upwind_dir_deg: upwind direction in degrees
"""
function create_kite_model(x, y, z, pos, upwind_dir_deg)
    kcu::KCU = KCU(se())
    s::KPS4 = KPS4(kcu)

    s.x = x
    s.y = y
    s.z = z

    KiteModels.set_v_wind_ground!(s, pos[begin+2]; upwind_dir=deg2rad(upwind_dir_deg))

    s.pos[end-2][begin] = pos[begin]
    s.pos[end-2][begin+1] = pos[begin+1]
    s.pos[end-2][begin+2] = pos[begin+2] 
    s
end  

"""
    obtain_results(s)

Return the tuple (roll, pitch, yaw, azimuth_north, elevation, heading) of the kite model s
"""
function obtain_results(s)
    roll, pitch, yaw = orient_euler(s)
    elevation = calc_elevation(s)
    azimuth_north = calc_azimuth_north(s)
    heading = calc_heading(s)
    return rad2deg(roll), rad2deg(pitch), rad2deg(yaw), rad2deg(azimuth_north), rad2deg(elevation), rad2deg(heading)
end

@testset verbose=true "Test roll, pitch, yaw, azimuth, elevation and heading..." begin

# Kite at an elevation of 45 degrees and with 0 roll pitch yaw and azimuth. Heading is then 180 degrees"
@testset "elevation 45" begin
    s = create_kite_model((0, 1, 0), (1, 0, 0), (0, 0, -1), # Orientation
                          (0, sqrt(2)/2, sqrt(2)/2))        # Pos ENU
    roll, pitch, yaw, azimuth_north, elevation, heading = obtain_results(s)

    @test isapprox(roll,            0,      atol=1e-4, rtol=1e-4)
    @test isapprox(pitch,           0,      atol=1e-4, rtol=1e-4)
    @test isapprox(yaw,             0,      atol=1e-4, rtol=1e-4)
    @test isapprox(azimuth_north,   0,      atol=1e-4, rtol=1e-4)
    @test isapprox(elevation,       45,     atol=1e-4, rtol=1e-4)
    @test isapprox(heading,         180,    atol=1e-4, rtol=1e-4)
end

# Kite at an elevation of 60 degrees and with 0 roll pitch yaw and azimuth. Heading is then 180 degrees"
@testset "elevation 60" begin
    s = create_kite_model((0, 1, 0), (1, 0, 0), (0, 0, -1), # Orientation
                          (0, 0.5, sqrt(3)/2))              # Pos ENU
    roll, pitch, yaw, azimuth_north, elevation, heading = obtain_results(s)

    @test isapprox(roll,            0,      atol=1e-4, rtol=1e-4)
    @test isapprox(pitch,           0,      atol=1e-4, rtol=1e-4)
    @test isapprox(yaw,             0,      atol=1e-4, rtol=1e-4)
    @test isapprox(azimuth_north,   0,      atol=1e-4, rtol=1e-4)
    @test isapprox(elevation,       60,     atol=1e-4, rtol=1e-4)
    @test isapprox(heading,         180,    atol=1e-4, rtol=1e-4)
end

# Kite at an elevation of -60 degrees and with 0 roll pitch yaw and azimuth. Heading is then 0 degrees"
@testset "elevation -60" begin
    s = create_kite_model((0, 1, 0), (1, 0, 0), (0, 0, -1), # Orientation
                          (0, 0.5, -sqrt(3)/2))             # Pos
    roll, pitch, yaw, azimuth_north, elevation, heading = obtain_results(s)

    @test isapprox(roll,            0,      atol=1e-4, rtol=1e-4)
    @test isapprox(pitch,           0,      atol=1e-4, rtol=1e-4)
    @test isapprox(yaw,             0,      atol=1e-4, rtol=1e-4)
    @test isapprox(azimuth_north,   0,      atol=1e-4, rtol=1e-4)
    @test isapprox(elevation,       -60,    atol=1e-4, rtol=1e-4)
    @test isapprox(heading,         0,      atol=1e-4, rtol=1e-4)
end

# Kite at an azimuth of 45 degrees and with 0 roll pitch yaw and elevation. Heading is then 270 degrees"
@testset "azimuth 45" begin
    s = create_kite_model((0, 1, 0), (1, 0, 0), (0, 0, -1), # Orientation
                          (-sqrt(2)/2, sqrt(2)/2, 0))       # Pos
    roll, pitch, yaw, azimuth_north, elevation, heading = obtain_results(s)

    @test isapprox(roll,            0,      atol=1e-4, rtol=1e-4)
    @test isapprox(pitch,           0,      atol=1e-4, rtol=1e-4)
    @test isapprox(yaw,             0,      atol=1e-4, rtol=1e-4)
    @test isapprox(azimuth_north,   45,     atol=1e-4, rtol=1e-4)
    @test isapprox(elevation,       0,      atol=1e-4, rtol=1e-4)
    @test isapprox(heading,         270,    atol=1e-4, rtol=1e-4)
end

# Kite at an azimuth of 60 degrees and with 0 roll pitch yaw and elevation. Heading is then 270 degrees"
@testset "azimuth 60" begin
    s = create_kite_model((0, 1, 0), (1, 0, 0), (0, 0, -1), # Orientation
                          (-sqrt(3)/2, 1/2, 0))             # Pos
    roll, pitch, yaw, azimuth_north, elevation, heading = obtain_results(s)

    @test isapprox(roll,            0,      atol=1e-4, rtol=1e-4)
    @test isapprox(pitch,           0,      atol=1e-4, rtol=1e-4)
    @test isapprox(yaw,             0,      atol=1e-4, rtol=1e-4)
    @test isapprox(azimuth_north,   60,     atol=1e-4, rtol=1e-4)
    @test isapprox(elevation,       0,      atol=1e-4, rtol=1e-4)
    @test isapprox(heading,         270,    atol=1e-4, rtol=1e-4)
end


# Kite at an azimuth of -60 degrees and with 0 roll pitch yaw and elevation. Heading is then 90 degrees"
@testset "azimuth -60" begin
    s = create_kite_model((0, 1, 0), (1, 0, 0), (0, 0, -1), # Orientation
                           (sqrt(3)/2, 1/2, 0))             # Pos
    roll, pitch, yaw, azimuth_north, elevation, heading = obtain_results(s)

    @test isapprox(roll,            0,      atol=1e-4, rtol=1e-4)
    @test isapprox(pitch,           0,      atol=1e-4, rtol=1e-4)
    @test isapprox(yaw,             0,      atol=1e-4, rtol=1e-4)
    @test isapprox(azimuth_north,   -60,    atol=1e-4, rtol=1e-4)
    @test isapprox(elevation,       0,      atol=1e-4, rtol=1e-4)
    @test isapprox(heading,         90,     atol=1e-4, rtol=1e-4)
end

# Kite at an azimuth of 60 degrees and with 45 degrees pitch, 0 roll yaw and elevation. Heading is then 319 degrees"
@testset "pitch azimuth" begin
    s = create_kite_model((0, sqrt(2)/2, sqrt(2)/2), (1, 0, 0), (0, sqrt(2)/2, -sqrt(2)/2),   # Orientation
                          (-sqrt(3)/2, 1/2, 0))                                                 # Pos
    roll, pitch, yaw, azimuth_north, elevation, heading = obtain_results(s)

    @test isapprox(roll,            0,      atol=1e-4, rtol=1e-4)
    @test isapprox(pitch,           45,     atol=1e-4, rtol=1e-4)
    @test isapprox(yaw,             0,      atol=1e-4, rtol=1e-4)
    @test isapprox(azimuth_north,   60,     atol=1e-4, rtol=1e-4)
    @test isapprox(elevation,       0,      atol=1e-4, rtol=1e-4)
    @test isapprox(heading,         319.1066053508691,    atol=1e-4, rtol=1e-4)
end

# Kite at an elevation of 60 degrees and with 45 degrees yaw, 0 roll pitch and elevation. Heading is then 229 degrees"
@testset "yaw elevation" begin
    s = create_kite_model((sqrt(2)/2, sqrt(2)/2, 0), (sqrt(2)/2, -sqrt(2)/2, 0), (0, 0, -1),    # Orientation
                          (0, 1/2, sqrt(3)/2))                                                  # Pos
    roll, pitch, yaw, azimuth_north, elevation, heading = obtain_results(s)

    @test isapprox(roll,            0,      atol=1e-4, rtol=1e-4)
    @test isapprox(pitch,           0,      atol=1e-4, rtol=1e-4)
    @test isapprox(yaw,             45,     atol=1e-4, rtol=1e-4)
    @test isapprox(azimuth_north,   0,      atol=1e-4, rtol=1e-4)
    @test isapprox(elevation,       60,     atol=1e-4, rtol=1e-4)
    @test isapprox(heading,         229.10660535086907,    atol=1e-4, rtol=1e-4)
end

# Kite at an elevation and azimuth of 45 degrees and with 0 roll pitch and yaw. Heading is then 234 degrees"
@testset "azimuth 45 elevation 45" begin
    s = create_kite_model((0, 1, 0), (1, 0, 0), (0, 0, -1), # Orientation
                          (-1/2, 1/2, sqrt(2)/2))           # Pos
    roll, pitch, yaw, azimuth_north, elevation, heading = obtain_results(s)

    @test isapprox(roll,            0,      atol=1e-4, rtol=1e-4)
    @test isapprox(pitch,           0,      atol=1e-4, rtol=1e-4)
    @test isapprox(yaw,             0,      atol=1e-4, rtol=1e-4)
    @test isapprox(azimuth_north,   45,     atol=1e-4, rtol=1e-4)
    @test isapprox(elevation,       45,     atol=1e-4, rtol=1e-4)
    @test isapprox(heading,         234.73561031724535,    atol=1e-4, rtol=1e-4)
end

# Kite at an elevation and azimuth of 45 degrees and 90 degrees roll. With 0 pitch and yaw. Heading is then 234 degrees and should not change due to roll"
@testset "roll azimuth elevation" begin
    s = create_kite_model((0, 1, 0), (0, 0, -1), (-1, 0, 0), # Orientation
                          (-1/2, 1/2, sqrt(2)/2))           # Pos
    roll, pitch, yaw, azimuth_north, elevation, heading = obtain_results(s)

    @test isapprox(roll,            90,     atol=1e-4, rtol=1e-4)
    @test isapprox(pitch,           0,      atol=1e-4, rtol=1e-4)
    @test isapprox(yaw,             0,      atol=1e-4, rtol=1e-4)
    @test isapprox(azimuth_north,   45,     atol=1e-4, rtol=1e-4)
    @test isapprox(elevation,       45,     atol=1e-4, rtol=1e-4)
    @test isapprox(heading,         234.73561031724535,    atol=1e-4, rtol=1e-4)
end

# Kite pointing towards zenith and z-axis towards groundstation at elevation and azimuth equal to 0, should give a heading of 0."
@testset "base orientation" begin
    s = create_kite_model((0, 0, 1), (-1, 0, 0), (0, -1, 0), # Orientation
                          (0, 1, 0))                         # Pos
    roll, pitch, yaw, azimuth_north, elevation, heading = obtain_results(s)

    @test isapprox(azimuth_north,   0,      atol=1e-4, rtol=1e-4)
    @test isapprox(elevation,       0,      atol=1e-4, rtol=1e-4)
    if heading > 359
        heading -= 360
    end
    @test isapprox(heading,         0,      atol=1e-4, rtol=1e-4)
end

# Kite pointing towards west and z-axis towards groundstation at elevation and azimuth equal to 0, should give a heading of 90."
@testset "base orientation heading 90" begin
    s = create_kite_model((-1, 0, 0), (0, 0, -1), (0, -1, 0), # Orientation
                          (0, 1, 0))                          # Pos
    roll, pitch, yaw, azimuth_north, elevation, heading = obtain_results(s)

    @test isapprox(azimuth_north,   0,      atol=1e-4, rtol=1e-4)
    @test isapprox(elevation,       0,      atol=1e-4, rtol=1e-4)
    @test isapprox(heading,         90,      atol=1e-4, rtol=1e-4)
end

# Kite is same place and orientation as base orientation, rotate the windframe x axis 45 degrees to the west"
@testset "upwind_dir dir 45" begin
    s = create_kite_model((0, 0, 1), (-1, 0, 0), (0, -1, 0), # Orientation
                          (0, 1, 0),                         # Pos
                          45+180)                            # upwind_dir  
    roll, pitch, yaw, azimuth_north, elevation, heading = obtain_results(s)

    @test isapprox(azimuth_north,   0,      atol=1e-4, rtol=1e-4)
    @test isapprox(elevation,       0,      atol=1e-4, rtol=1e-4)
    @test isapprox(heading,         0,      atol=1e-4, rtol=1e-4) || isapprox(heading, 360, atol=1e-4, rtol=1e-4)
end

# Kite is same place and orientation as base orientation, rotate the windframe x axis 60 degrees to the west"
@testset "upwind_dir dir 60" begin
    s = create_kite_model((0, 0, 1), (-1, 0, 0), (0, -1, 0), # Orientation
                          (0, 1, 0),                         # Pos
                          60+180)                            # upwind_dir  
    roll, pitch, yaw, azimuth_north, elevation, heading = obtain_results(s)

    @test isapprox(azimuth_north,   0,      atol=1e-4, rtol=1e-4)
    @test isapprox(elevation,       0,      atol=1e-4, rtol=1e-4)
    @test isapprox(heading,         0,      atol=1e-4, rtol=1e-4)
end

# Kite is same place and orientation as base orientation, rotate the windframe x axis 60 degrees to the east"
@testset "upwind_dir dir -60" begin
    s = create_kite_model((0, 0, 1), (-1, 0, 0), (0, -1, 0), # Orientation
                          (0, 1, 0),                         # Pos
                          -60+180)                           # upwind_dir  
    roll, pitch, yaw, azimuth_north, elevation, heading = obtain_results(s)
    
    @test isapprox(azimuth_north,   0,      atol=1e-4, rtol=1e-4)
    @test isapprox(elevation,       0,      atol=1e-4, rtol=1e-4)
    @test isapprox(heading,         0,      atol=1e-4, rtol=1e-4)
end

# Kite is same place and orientation as base orientation, rotate the windframe x axis 45 degrees to the west"
@testset "upwind_dir dir 45 azimuth wind" begin
    s = create_kite_model((0, 0, 1), (-1, 0, 0), (0, -1, 0), # Orientation
                          (0, 1, 0),                         # Pos
                          45+180)                            # upwind_dir  
    azimuth = rad2deg(calc_azimuth(s))
    azimuth_north = rad2deg(calc_azimuth_north(s))
    @test isapprox(azimuth,         45,      atol=1e-4, rtol=1e-4)
    @test isapprox(azimuth_north,   0,      atol=1e-4, rtol=1e-4)
end

# Kite is same place and orientation as base orientation, rotate the windframe x axis 60 degrees to the west"
@testset "upwind_dir dir 60 azimuth wind" begin
    s = create_kite_model((0, 0, 1), (-1, 0, 0), (0, -1, 0), # Orientation
                          (0, 1, 0),                         # Pos
                          60+180)                            # upwind_dir  
    upwind_dir_ = rad2deg(upwind_dir(s))
    azimuth = rad2deg(calc_azimuth(s))
    azimuth_north = rad2deg(calc_azimuth_north(s))

    @test isapprox(azimuth,         60,      atol=1e-4, rtol=1e-4)
    @test isapprox(azimuth_north,   0,      atol=1e-4, rtol=1e-4)
end

# Kite is same place and orientation as base orientation, rotate the windframe x axis 60 degrees to the east"
@testset "upwind_dir dir -60 azimuth wind" begin
    s = create_kite_model((0, 0, 1), (-1, 0, 0), (0, -1, 0), # Orientation
                          (0, 1, 0),                         # Pos
                          -60+180)                           # upwind_dir  
    azimuth = rad2deg(calc_azimuth(s))
    azimuth_north = rad2deg(calc_azimuth_north(s))

    @test isapprox(azimuth,         -60,      atol=1e-4, rtol=1e-4)
    @test isapprox(azimuth_north,   0,      atol=1e-4, rtol=1e-4)

end
end
nothing
    
calc_azimuth