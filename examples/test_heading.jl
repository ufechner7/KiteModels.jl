using KiteUtils
using KiteModels
using KitePodModels

using LinearAlgebra
using StaticArrays
using Test

 """
     fromW2SE(vector, elevation, azimuth)

 Transform a (velocity-) vector (x,y,z) from Wind to Small Earth reference frame .
 """
function fromW2SE2(vector, elevation, azimuth)
    rotate_first_step = @SMatrix[0  0  1;
                                 0  1  0;
                                -1  0  0]
    rotate_elevation = @SMatrix[cos(elevation) 0 sin(elevation);
                                0              1         0;
                             -sin(elevation)   0   cos(elevation)]
    rotate_azimuth = @SMatrix[1         0       0;
                              0  cos(-azimuth)   -sin(-azimuth);
                              0  sin(-azimuth)    cos(-azimuth)]
    rotate_elevation * rotate_azimuth * rotate_first_step * vector
end

function calc_heading_w2(orientation, down_wind_direction = pi/2.0)
    # create a unit heading vector in the xsense reference frame
    heading_sensor =  SVector(1, 0, 0)
    # rotate headingSensor to the Earth Xsens reference frame
    headingEX = fromKS2EX(heading_sensor, orientation)
    # rotate headingEX to earth groundstation reference frame
    headingEG = fromEX2EG(headingEX)
    # rotate headingEG to headingW and convert to 2d HeadingW vector
    fromEG2W(headingEG, down_wind_direction)
end

function calc_heading2(orientation, elevation, azimuth; upwind_dir=-pi/2, respos=true)
    down_wind_direction = wrap2pi(upwind_dir + π)
    headingSE = fromW2SE2(calc_heading_w2(orientation, down_wind_direction), elevation, azimuth)
    angle = atan(headingSE.y, headingSE.x)
    if angle < 0 && respos
        angle += 2π
    end
    angle
end

function calc_heading2(s::KiteModels.AKM; upwind_dir_=upwind_dir(s))
    orientation = orient_euler(s)
    elevation = calc_elevation(s)
    # use azimuth in wind reference frame
    azimuth = calc_azimuth(s)
    calc_heading2(orientation, elevation, azimuth; upwind_dir=upwind_dir_)
end


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

    KiteModels.set_v_wind_ground!(s, pos[begin+2], deg2rad(upwind_dir_deg))

    s.pos[end-2][begin] = pos[begin]
    s.pos[end-2][begin+1] = pos[begin+1]
    s.pos[end-2][begin+2] = pos[begin+2] 
    s
end  

@testset verbose=true "Test heading..." begin

# Kite at an elevation of 45 degrees and with 0 roll pitch yaw and azimuth. Heading is then 180 degrees"
@testset "elevation 45" begin
    s = create_kite_model((0, 1, 0), (1, 0, 0), (0, 0, -1), # Orientation
                          (0, sqrt(2)/2, sqrt(2)/2))        # Pos ENU
    heading = rad2deg(calc_heading2(s))
    @test isapprox(heading,         180,    atol=1e-4, rtol=1e-4)
end
@testset "elevation 60" begin
    s = create_kite_model((0, 1, 0), (1, 0, 0), (0, 0, -1), # Orientation
                          (0, 0.5, sqrt(3)/2))              # Pos ENU
    heading = rad2deg(calc_heading2(s))
    @test isapprox(heading,         180,    atol=1e-4, rtol=1e-4)
end
# Kite at an elevation of -60 degrees and with 0 roll pitch yaw and azimuth. Heading is then 0 degrees"
@testset "elevation -60" begin
    s = create_kite_model((0, 1, 0), (1, 0, 0), (0, 0, -1), # Orientation
                          (0, 0.5, -sqrt(3)/2))             # Pos
    heading = rad2deg(calc_heading2(s))
    @test isapprox(heading,         0,      atol=1e-4, rtol=1e-4)
end
# Kite at an azimuth of 45 degrees and with 0 roll pitch yaw and elevation. Heading is then 270 degrees"
@testset "azimuth 45" begin
    s = create_kite_model((0, 1, 0), (1, 0, 0), (0, 0, -1), # Orientation
                          (-sqrt(2)/2, sqrt(2)/2, 0))       # Pos
    heading = rad2deg(calc_heading2(s))
    @test isapprox(heading,         270,    atol=1e-4, rtol=1e-4)
end
# Kite at an azimuth of 60 degrees and with 0 roll pitch yaw and elevation. Heading is then 270 degrees"
@testset "azimuth 60" begin
    s = create_kite_model((0, 1, 0), (1, 0, 0), (0, 0, -1), # Orientation
                          (-sqrt(3)/2, 1/2, 0))             # Pos
    heading = rad2deg(calc_heading2(s))
    @test isapprox(heading,         270,    atol=1e-4, rtol=1e-4)
end
# Kite at an azimuth of -60 degrees and with 0 roll pitch yaw and elevation. Heading is then 90 degrees"
@testset "azimuth_north -60" begin
    s = create_kite_model((0, 1, 0), (1, 0, 0), (0, 0, -1), # Orientation
                           (sqrt(3)/2, 1/2, 0))             # Pos
    heading = rad2deg(calc_heading2(s))
    @test isapprox(heading,         90,     atol=1e-4, rtol=1e-4)
end
# Kite at an azimuth of -60 degrees and with 0 roll pitch yaw and elevation. Heading is then 90 degrees"
@testset "azimuth_north -60" begin
    s = create_kite_model((0, 1, 0), (1, 0, 0), (0, 0, -1), # Orientation
                           (sqrt(3)/2, 1/2, 0))             # Pos
    heading = rad2deg(calc_heading2(s))
    @test isapprox(heading,         90,     atol=1e-4, rtol=1e-4)
end
# Kite at an azimuth of 60 degrees and with 45 degrees pitch, 0 roll yaw and elevation. Heading is then 315 degrees"
@testset "pitch 45 azimuth_north 60" begin
    s = create_kite_model((0, sqrt(2)/2, sqrt(2)/2), (1, 0, 0), (0, sqrt(2)/2, -sqrt(2)/2),   # Orientation
                          (-sqrt(3)/2, 1/2, 0))                                               # Pos
    heading = rad2deg(calc_heading2(s))
    println(heading)
    @test_broken isapprox(heading,         315,    atol=1e-4, rtol=1e-4)
end
end
nothing