using Pkg

using KiteUtils
using KiteModels
using KitePodModels

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
    return rad2deg(roll), rad2deg(pitch), rad2deg(yaw), rad2deg(azimuth), rad2deg(elevation), rad2deg(heading)
end

function print_results(test, roll, pitch, yaw, elevation, azimuth, heading)
    print(test, "\n")
    @printf("Roll: %0.2f, Pitch: %0.2f, Yaw: %0.2f\n", roll, pitch, yaw)
    @printf("Elevation: %0.2f, Azimuth: %0.2f, Heading: %0.2f\n", elevation, azimuth, heading)
    
end

# Kite at an elevation of 45 degrees and with 0 roll pitch yaw and azimuth. Heading is then 180 degrees
@testset "elevation 45" begin
    s = create_kite_model((0, 1, 0), (1, 0, 0), (0, 0, -1), # Orientation
                          (0, sqrt(2)/2, sqrt(2)/2))        # Pos
    roll, pitch, yaw, azimuth, elevation, heading = obtain_results(s)

    print_results("test elevation 45", roll, pitch, yaw, elevation, azimuth, heading)

    @test isapprox(roll,        0,      atol=1e-4, rtol=1e-4)
    @test isapprox(pitch,       0,      atol=1e-4, rtol=1e-4)
    @test isapprox(yaw,         0,      atol=1e-4, rtol=1e-4)
    @test isapprox(azimuth,     0,      atol=1e-4, rtol=1e-4)
    @test isapprox(elevation,   45,     atol=1e-4, rtol=1e-4)
    @test isapprox(heading,     180,    atol=1e-4, rtol=1e-4)
end

# Kite at an elevation of 60 degrees and with 0 roll pitch yaw and azimuth. Heading is then 180 degrees
@testset "elevation 60" begin
    s = create_kite_model((0, 1, 0), (1, 0, 0), (0, 0, -1), # Orientation
                          (0, 0.5, sqrt(3)/2))              # Pos
    roll, pitch, yaw, azimuth, elevation, heading = obtain_results(s)
    
    print_results("test elevation 60", roll, pitch, yaw, elevation, azimuth, heading)

    @test isapprox(roll,        0,      atol=1e-4, rtol=1e-4)
    @test isapprox(pitch,       0,      atol=1e-4, rtol=1e-4)
    @test isapprox(yaw,         0,      atol=1e-4, rtol=1e-4)
    @test isapprox(azimuth,     0,      atol=1e-4, rtol=1e-4)
    @test isapprox(elevation,   60,     atol=1e-4, rtol=1e-4)
    @test isapprox(heading,     180,    atol=1e-4, rtol=1e-4)
end

# Kite at an elevation of -60 degrees and with 0 roll pitch yaw and azimuth. Heading is then 0 degrees"
@testset "elevation -60" begin
    s = create_kite_model((0, 1, 0), (1, 0, 0), (0, 0, -1), # Orientation
                          (0, 0.5, -sqrt(3)/2))             # Pos
    roll, pitch, yaw, azimuth, elevation, heading = obtain_results(s)
    
    print_results("test elevation -60", roll, pitch, yaw, elevation, azimuth, heading)

    @test isapprox(roll,        0,      atol=1e-4, rtol=1e-4)
    @test isapprox(pitch,       0,      atol=1e-4, rtol=1e-4)
    @test isapprox(yaw,         0,      atol=1e-4, rtol=1e-4)
    @test isapprox(azimuth,     0,      atol=1e-4, rtol=1e-4)
    @test isapprox(elevation,   -60,    atol=1e-4, rtol=1e-4)
    @test isapprox(heading,     0,      atol=1e-4, rtol=1e-4)
end

# Kite at an azimuth of 45 degrees and with 0 roll pitch yaw and elevation. Heading is then 270 degrees"
@testset "azimuth 45" begin
    s = create_kite_model((0, 1, 0), (1, 0, 0), (0, 0, -1), # Orientation
                          (-sqrt(2)/2, sqrt(2)/2, 0))       # Pos
    roll, pitch, yaw, azimuth, elevation, heading = obtain_results(s)
    
    print_results("test azimuth 45", roll, pitch, yaw, elevation, azimuth, heading)

    @test isapprox(roll,        0,      atol=1e-4, rtol=1e-4)
    @test isapprox(pitch,       0,      atol=1e-4, rtol=1e-4)
    @test isapprox(yaw,         0,      atol=1e-4, rtol=1e-4)
    @test isapprox(azimuth,     45,     atol=1e-4, rtol=1e-4)
    @test isapprox(elevation,   0,      atol=1e-4, rtol=1e-4)
    @test isapprox(heading,     270,    atol=1e-4, rtol=1e-4)
end

# Kite at an azimuth of 60 degrees and with 0 roll pitch yaw and elevation. Heading is then 270 degrees"
@testset "azimuth 60" begin
    s = create_kite_model((0, 1, 0), (1, 0, 0), (0, 0, -1), # Orientation
                          (-sqrt(3)/2, 1/2, 0))             # Pos
    roll, pitch, yaw, azimuth, elevation, heading = obtain_results(s)
    
    print_results("test azimuth 60", roll, pitch, yaw, elevation, azimuth, heading)

    @test isapprox(roll,        0,      atol=1e-4, rtol=1e-4)
    @test isapprox(pitch,       0,      atol=1e-4, rtol=1e-4)
    @test isapprox(yaw,         0,      atol=1e-4, rtol=1e-4)
    @test isapprox(azimuth,     60,     atol=1e-4, rtol=1e-4)
    @test isapprox(elevation,   0,      atol=1e-4, rtol=1e-4)
    @test isapprox(heading,     270,    atol=1e-4, rtol=1e-4)
end


# Kite at an azimuth of -60 degrees and with 0 roll pitch yaw and elevation. Heading is then 90 degrees"
@testset "azimuth -60" begin
    s = create_kite_model((0, 1, 0), (1, 0, 0), (0, 0, -1), # Orientation
                           (sqrt(3)/2, 1/2, 0))             # Pos
    roll, pitch, yaw, azimuth, elevation, heading = obtain_results(s)
    
    print_results("test azimuth -60", roll, pitch, yaw, elevation, azimuth, heading)

    @test isapprox(roll,        0,      atol=1e-4, rtol=1e-4)
    @test isapprox(pitch,       0,      atol=1e-4, rtol=1e-4)
    @test isapprox(yaw,         0,      atol=1e-4, rtol=1e-4)
    @test isapprox(azimuth,     -60,    atol=1e-4, rtol=1e-4)
    @test isapprox(elevation,   0,      atol=1e-4, rtol=1e-4)
    @test isapprox(heading,     90,    atol=1e-4, rtol=1e-4)
end

# Kite at an azimuth of 60 degrees and with 45 degrees pitch, 0 roll yaw and elevation. Heading is then 315 degrees"
@testset "pitch azimuth" begin
    s = create_kite_model((0, sqrt(2)/2, sqrt(2)/2), (1, 0, 0), (0, sqrt(2)/2, -sqrt(2)/2),   # Orientation
                          (-sqrt(3)/2, 1/2, 0))                                                 # Pos
    roll, pitch, yaw, azimuth, elevation, heading = obtain_results(s)
    
    print_results("test ptich azimuth", roll, pitch, yaw, elevation, azimuth, heading)

    @test isapprox(roll,        0,      atol=1e-4, rtol=1e-4)
    @test isapprox(pitch,       45,     atol=1e-4, rtol=1e-4)
    @test isapprox(yaw,         0,      atol=1e-4, rtol=1e-4)
    @test isapprox(azimuth,     60,     atol=1e-4, rtol=1e-4)
    @test isapprox(elevation,   0,      atol=1e-4, rtol=1e-4)
    @test isapprox(heading,     315,    atol=1e-4, rtol=1e-4)
end

# Kite at an elevation of 60 degrees and with 45 degrees yaw, 0 roll pitch and elevation. Heading is then 225 degrees"
@testset "yaw elevation" begin
    s = create_kite_model((sqrt(2)/2, sqrt(2)/2, 0), (sqrt(2)/2, sqrt(2)/2, 0), (0, 0, -1),    # Orientation
                          (0, 1/2, sqrt(3)/2))                                                  # Pos
    roll, pitch, yaw, azimuth, elevation, heading = obtain_results(s)
    
    print_results("test yaw elevation", roll, pitch, yaw, elevation, azimuth, heading)

    @test isapprox(roll,        0,      atol=1e-4, rtol=1e-4)
    @test isapprox(pitch,       0,      atol=1e-4, rtol=1e-4)
    @test isapprox(yaw,         45,     atol=1e-4, rtol=1e-4)
    @test isapprox(azimuth,     0,      atol=1e-4, rtol=1e-4)
    @test isapprox(elevation,   60,     atol=1e-4, rtol=1e-4)
    @test isapprox(heading,     135,    atol=1e-4, rtol=1e-4)
end

# Kite at an elevation and azimuth of 45 degrees and with 0 roll pitch and yaw. Heading is then 225 degrees"
@testset "azimuth 45 elevation 45" begin
    s = create_kite_model((0, 1, 0), (1, 0, 0), (0, 0, -1), # Orientation
                          (-1/2, 1/2, sqrt(2)/2))           # Pos
    roll, pitch, yaw, azimuth, elevation, heading = obtain_results(s)
    
    print_results("test azimuth 45 elevation 45", roll, pitch, yaw, elevation, azimuth, heading)

    @test isapprox(roll,        0,      atol=1e-4, rtol=1e-4)
    @test isapprox(pitch,       0,      atol=1e-4, rtol=1e-4)
    @test isapprox(yaw,         0,      atol=1e-4, rtol=1e-4)
    @test isapprox(azimuth,     45,     atol=1e-4, rtol=1e-4)
    @test isapprox(elevation,   45,     atol=1e-4, rtol=1e-4)
    @test isapprox(heading,     225,    atol=1e-4, rtol=1e-4)
end

# Kite at an elevation and azimuth of 45 degrees and 90 degrees roll. With 0 pitch and yaw. Heading is then 225 degrees and should not change due to roll"
@testset "roll azimuth elevation" begin
    s = create_kite_model((0, 1, 0), (0, 0, -1), (-1, 0, 0), # Orientation
                          (-1/2, 1/2, sqrt(2)/2))           # Pos
    roll, pitch, yaw, azimuth, elevation, heading = obtain_results(s)
    
    print_results("test azimuth elevation", roll, pitch, yaw, elevation, azimuth, heading)

    @test isapprox(roll,        90,     atol=1e-4, rtol=1e-4)
    @test isapprox(pitch,       0,      atol=1e-4, rtol=1e-4)
    @test isapprox(yaw,         0,      atol=1e-4, rtol=1e-4)
    @test isapprox(azimuth,     45,     atol=1e-4, rtol=1e-4)
    @test isapprox(elevation,   45,     atol=1e-4, rtol=1e-4)
    @test isapprox(heading,     225,    atol=1e-4, rtol=1e-4)
end

# Kite pointing towards zenith and z-axis towards groundstation at elevation and azimuth equal to 0, should give a heading of 0."
@testset "base orientation" begin
    s = create_kite_model((0, 0, 1), (-1, 0, 0), (0, -1, 0), # Orientation
                          (0, 1, 0))                         # Pos
    roll, pitch, yaw, azimuth, elevation, heading = obtain_results(s)
    
    print_results("test base orientation", roll, pitch, yaw, elevation, azimuth, heading)

    @test isapprox(azimuth,     0,      atol=1e-4, rtol=1e-4)
    @test isapprox(elevation,   0,      atol=1e-4, rtol=1e-4)
    @test isapprox(heading,     0,      atol=1e-4, rtol=1e-4)
end

# Kite pointing towards west and z-axis towards groundstation at elevation and azimuth equal to 0, should give a heading of 90."
@testset "base orientation heading 90" begin
    s = create_kite_model((-1, 0, 0), (0, 0, -1), (0, -1, 0), # Orientation
                          (0, 1, 0))                          # Pos
    roll, pitch, yaw, azimuth, elevation, heading = obtain_results(s)
    
    print_results("test orientation heading 90", roll, pitch, yaw, elevation, azimuth, heading)

    @test isapprox(azimuth,     0,      atol=1e-4, rtol=1e-4)
    @test isapprox(elevation,   0,      atol=1e-4, rtol=1e-4)
    @test isapprox(heading,     0,      atol=1e-4, rtol=1e-4)
end

# Kite is same place and orientation as base orientation, rotate the windframe x axis 45 degrees to the west"
@testset "wind dir 45" begin
    s = create_kite_model((0, 0, 1), (-1, 0, 0), (0, -1, 0), # Orientation
                          (0, 1, 0),                         # Pos
                          45+180)                                # Wind dir  
    roll, pitch, yaw, azimuth, elevation, heading = obtain_results(s)
    
    print_results("test wind 45", roll, pitch, yaw, elevation, azimuth, heading)

    @test isapprox(azimuth,     -45,    atol=1e-4, rtol=1e-4)
    @test isapprox(elevation,   0,      atol=1e-4, rtol=1e-4)
    @test isapprox(heading,     0,      atol=1e-4, rtol=1e-4)
end

# Kite is same place and orientation as base orientation, rotate the windframe x axis 60 degrees to the west"
@testset "wind dir 60" begin
    s = create_kite_model((0, 0, 1), (-1, 0, 0), (0, -1, 0), # Orientation
                          (0, 1, 0),                         # Pos
                          60+180)                                 # Wind dir  
    roll, pitch, yaw, azimuth, elevation, heading = obtain_results(s)
    
    print_results("test wind 60", roll, pitch, yaw, elevation, azimuth, heading)

    @test isapprox(azimuth,     -60,    atol=1e-4, rtol=1e-4)
    @test isapprox(elevation,   0,      atol=1e-4, rtol=1e-4)
    @test isapprox(heading,     0,      atol=1e-4, rtol=1e-4)
end

# Kite is same place and orientation as base orientation, rotate the windframe x axis 60 degrees to the east"
@testset "wind dir -60" begin
    s = create_kite_model((0, 0, 1), (-1, 0, 0), (0, -1, 0), # Orientation
                          (0, 1, 0),                         # Pos
                          -60+180)                           # Wind dir  
    roll, pitch, yaw, azimuth, elevation, heading = obtain_results(s)
    
    print_results("test wind -60", roll, pitch, yaw, elevation, azimuth, heading)

    @test isapprox(azimuth,     60,      atol=1e-4, rtol=1e-4)
    @test isapprox(elevation,   0,     atol=1e-4, rtol=1e-4)
    @test isapprox(heading,     0,    atol=1e-4, rtol=1e-4)
end
    