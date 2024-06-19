#= MIT License

Copyright (c) 2020, 2021, 2022 Uwe Fechner

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE. =#

#= Model of a kite-power system in implicit form: residual = f(y, yd)

This model implements a 3D mass-spring system with reel-out. It uses six tether segments (the number can be
configured in the file data/settings.yaml). The kite is modelled using 4 point masses and 3*n aerodynamic 
surfaces. The spring constant and the damping decrease with the segment length. The aerodynamic kite forces
are acting on three of the four kite point masses. 

Four point kite model, included from KiteModels.jl.

Scientific background: http://arxiv.org/abs/1406.6218 =#

# Array of connections of bridlepoints.
# First point, second point, unstressed length.
const SPRINGS_INPUT_3L = [1.      4.  -1. # s1: E, A
                        1.      2.  -1. # s2, E, C
                        1.      3.  -1. # s3, E, D
                        2.      3.  -1. # s4, C, D
                        2.      4.  -1. # s5, C, A
                        3.      4.  -1. # s6, D, A
                        ]
# E = 1, C = 2, D = 3, A = 4
# E = segments*3+1, C = segments*3+2, D = segments*3+3, A = segments*3+4


const KITE_SPRINGS_3L = 6
const KITE_PARTICLES_3L = 4
const KITE_ANGLE_3L = 0.0
# struct, defining the phyical parameters of one spring
# @with_kw struct Spring{I, S}
#     p1::I = 1         # number of the first point
#     p2::I = 2         # number of the second point
#     length::S = 1.0   # current unstressed spring length
#     c_spring::S = 1.0 # spring constant [N/m]
#     damping::S  = 0.1 # damping coefficent [Ns/m]
# end

# const SP = Spring{Int16, SimFloat}
# const PRE_STRESS  = 0.9998   # Multiplier for the initial spring lengths.
# const KS = deg2rad(16.565 * 1.064 * 0.875 * 1.033 * 0.9757 * 1.083)  # max steering
# const DRAG_CORR = 0.93       # correction of the drag for the 4-point model
# function zero(::Type{SP})
#     SP(0,0,0,0,0)
# end

"""
    mutable struct KPS4_3L{S, T, P, Q, SP} <: AbstractKiteModel

State of the kite power system, using a 3 point kite model and three steering lines to the ground. Parameters:
- S: Scalar type, e.g. SimFloat
  In the documentation mentioned as Any, but when used in this module it is always SimFloat and not Any.
- T: Vector type, e.g. MVector{3, SimFloat}
- P: number of points of the system, segments+1
- Q: number of springs in the system, P-1
- SP: struct type, describing a spring
Normally a user of this package will not have to access any of the members of this type directly,
use the input and output functions instead.

$(TYPEDFIELDS)
"""
@with_kw mutable struct KPS4_3L{S, T, P, Q, SP} <: AbstractKiteModel
    "Reference to the settings struct"
    set::Settings = se()
    "Reference to the atmospheric model as implemented in the package AtmosphericModels"
    am::AtmosphericModel = AtmosphericModel()
    "Reference to the motor models as implemented in the package WinchModels. index 1: middle motor, index 2: left motor, index 3: right motor"
    motors::Vector{AbstractWinchModel} = [AsyncMachine() for _ in 1:3]
    "Iterations, number of calls to the function residual!"
    iter:: Int64 = 0
    "Function for calculation the lift coefficent, using a spline based on the provided value pairs."
    calc_cl = Spline1D(se().alpha_cl, se().cl_list)
    "Function for calculation the drag coefficent, using a spline based on the provided value pairs."
    calc_cd = Spline1D(se().alpha_cd, se().cd_list)   
    "wind vector at the height of the kite" 
    v_wind::T =           zeros(S, 3)
    "wind vector at reference height" 
    v_wind_gnd::T =       zeros(S, 3)
    "wind vector used for the calculation of the tether drag"
    v_wind_tether::T =    zeros(S, 3)
    "apparent wind vector at the kite"
    v_apparent::T =       zeros(S, 3)
    "drag force of kite and bridle; output of calc_aero_forces!"
    drag_force::T =       zeros(S, 3)
    "lift force of the kite; output of calc_aero_forces!"
    lift_force::T =       zeros(S, 3)    
    "spring force of the current tether segment, output of calc_particle_forces!"
    spring_force::T =     zeros(S, 3)
    "last winch force"
    last_force::T =       zeros(S, 3)
    "a copy of the residual one (pos,vel) for debugging and unit tests"    
    res1::SVector{P, KVec3} = zeros(SVector{P, KVec3})
    "a copy of the residual two (vel,acc) for debugging and unit tests"
    res2::SVector{P, KVec3} = zeros(SVector{P, KVec3})
    "a copy of the actual positions as output for the user"
    pos::SVector{P, KVec3} = zeros(SVector{P, KVec3})
    "velocity vector of the kite"
    vel_kite::T =          zeros(S, 3)
    "unstressed segment lengths of the three tethers [m]"
    segment_lengths::T =           zeros(S, 3)
    "lift coefficient of the kite, depending on the angle of attack"
    param_cl::S =         0.2
    "drag coefficient of the kite, depending on the angle of attack"
    param_cd::S =         1.0
    "azimuth angle in radian; inital value is zero"
    psi::S =              zero(S)
    # "depower angle [deg]"
    # alpha_depower::S =     0.0
    "relative start time of the current time interval"
    t_0::S =               0.0
    "reel out speed of the winch"
    reel_out_speeds::T =        zeros(S, 3)
    "reel out speed at the last time step"
    last_reel_out_speeds::T =   zeros(S, 3)
    "unstretched tether length"
    l_tethers::T =          zeros(S, 3)
    "lengths of the connections of the steering tethers to the kite"
    l_connections::T =      zeros(S, 2)
    "air density at the height of the kite"
    rho::S =               0.0
    # "actual relative depower setting,  must be between    0 .. 1.0"
    # depower::S =           0.0
    # "actual relative steering setting, must be between -1.0 .. 1.0"
    # steering::S =          0.0
    "multiplier for the stiffniss of tether and bridle"
    stiffness_factor::S =  1.0
    "initial masses of the point masses"
    initial_masses::MVector{P, S} = ones(P)
    "current masses, depending on the total tether length"
    masses::MVector{P, S}         = zeros(P)
    "vector of the springs, defined as struct"
    springs::MVector{Q, SP}       = zeros(SP, Q)
    "vector of the forces, acting on the particles"
    forces::SVector{P, KVec3} = zeros(SVector{P, KVec3})
    "synchronous speeds of the motors"
    sync_speeds::T =        zeros(S, 3)
    "x vector of kite reference frame"
    e_x::T =                 zeros(S, 3)
    "y vector of kite reference frame"
    e_y::T =                 zeros(S, 3)
    "z vector of kite reference frame"
    e_z::T =                 zeros(S, 3)
    "Point number of E"
    num_E::Int64 =           0
    "Point number of C"
    num_C::Int64 =           0
    "Point number of D"
    num_D::Int64 =           0
    "Point number of A"
    num_A::Int64 =           0
end

"""
    clear!(s::KPS4_3L)

Initialize the kite power model.
"""
function clear!(s::KPS4_3L)
    s.t_0 = 0.0                              # relative start time of the current time interval
    s.reel_out_speeds = zeros(3)
    s.last_reel_out_speeds = zeros(3)
    s.v_wind_gnd    .= [s.set.v_wind, 0, 0]    # wind vector at reference height
    s.v_wind_tether .= [s.set.v_wind, 0, 0]
    s.v_apparent    .= [s.set.v_wind, 0, 0]
    height = sin(deg2rad(s.set.elevation)) * (s.set.l_tether)
    s.v_wind .= s.v_wind_gnd * calc_wind_factor(s.am, height)

    s.l_tethers = [s.set.l_tether for _ in 1:3]
    s.segment_lengths = s.l_tethers ./ s.set.segments
    init_masses!(s)
    init_springs!(s)
    for i in 1:s.set.segments + KiteModels.KITE_PARTICLES + 1 
        s.forces[i] .= zeros(3)
    end
    s.drag_force .= [0.0, 0, 0]
    s.lift_force .= [0.0, 0, 0]
    s.rho = s.set.rho_0
    s.calc_cl = Spline1D(s.set.alpha_cl, s.set.cl_list)
    s.calc_cd = Spline1D(s.set.alpha_cd, s.set.cd_list) 
end

"""
    set_vs_reel_out!(s::KPS4_3L, reel_out_speeds, t_0, period_time = 1.0 / s.set.sample_freq)

Setter for the reel-out velocities of the three tethers. This manouvers the kite. Must be called on every timestep (before each simulation).
It also updates the tether lengths, therefore it must be called even if `reel_out_speeds` has
not changed.

- t_0 the start time of the next timestep relative to the start of the simulation [s]
"""
function set_reel_out_speeds!(s::KPS4_3L, reel_out_speeds::Vector{SimFloat}, t_0, period_time = 1.0 / s.set.sample_freq)
    s.sync_speeds .= reel_out_speeds
    s.last_reel_out_speeds = s.reel_out_speeds
    s.t_0 = t_0
end

function KPS4_3L()
    set = se()
    s = KPS4_3L{SimFloat, KVec3, set.segments*3+2+KITE_PARTICLES, set.segments*3+KITE_SPRINGS_3L, SP}()
    s.set = set
    # E = segments*3+1, C = segments*3+2, D = segments*3+3, A = segments*3+4
    s.num_E = s.set.segments*3+3
    s.num_C = s.set.segments*3+3+1
    s.num_D = s.set.segments*3+3+2
    s.num_A = s.set.segments*3+3+3
    s.calc_cl = Spline1D(s.set.alpha_cl, s.set.cl_list)
    s.calc_cd = Spline1D(s.set.alpha_cd, s.set.cd_list)       
    clear!(s)
    return s
end

""" 
    calc_particle_forces!(s::KPS4_3L, pos1, pos2, vel1, vel2, spring, segments, d_tether, rho, i)

Calculate the drag force and spring force of the tether segment, defined by the parameters pos1, pos2, vel1 and vel2
and distribute it equally on the two particles, that are attached to the segment.
The result is stored in the array s.forces. 
"""
@inline function calc_particle_forces!(s::KPS4_3L, pos1, pos2, vel1, vel2, spring, segments, d_tether, rho, i)
    l_0 = spring.length # Unstressed length
    k = spring.c_spring * s.stiffness_factor  # Spring constant
    c = spring.damping                        # Damping coefficient    
    segment = pos1 - pos2
    rel_vel = vel1 - vel2
    av_vel = 0.5 * (vel1 + vel2)
    norm1 = norm(segment)
    unit_vector = segment / norm1

    k1 = 0.25 * k # compression stiffness kite segments
    k2 = 0.1 * k  # compression stiffness tether segments
    c1 = 6.0 * c  # damping kite segments
    spring_vel   = unit_vector ⋅ rel_vel
    if (norm1 - l_0) > 0.0
        if i > segments  # kite springs
             s.spring_force .= (k *  (norm1 - l_0) + (c1 * spring_vel)) * unit_vector 
        else
             s.spring_force .= (k *  (norm1 - l_0) + (c * spring_vel)) * unit_vector
        end
    elseif i > segments # kite spring
        s.spring_force .= (k1 *  (norm1 - l_0) + (c * spring_vel)) * unit_vector
    else
        s.spring_force .= (k2 *  (norm1 - l_0) + (c * spring_vel)) * unit_vector
    end

    s.v_apparent .= s.v_wind_tether - av_vel
    if s.set.version == 1
        area = norm1 * d_tether
    else
        area = norm1 * s.set.d_line * 0.001
    end
    v_app_perp = s.v_apparent - s.v_apparent ⋅ unit_vector * unit_vector
    half_drag_force = (-0.25 * rho * s.set.cd_tether * norm(v_app_perp) * area) * v_app_perp 

    @inbounds s.forces[spring.p1] .+= half_drag_force + s.spring_force
    @inbounds s.forces[spring.p2] .+= half_drag_force - s.spring_force
    if i == 1 s.last_force .= s.forces[spring.p1] end
    nothing
end

function calc_kite_ref_frame!(s::KPS4_3L, E, C, D)
    P_c = 0.5 .* (C+D)
    s.e_y = normalize(C - D)
    s.e_z = normalize(E - P_c)
    s.e_x = cross(e_y, e_z)
    return e_x, e_y, e_z
end

"""
    calc_aero_forces!(s::KPS4_3L, pos, vel, rho, alpha_depower, rel_steering)

Calculates the aerodynamic forces acting on the kite particles.

Parameters:
- pos:              vector of the particle positions
- vel:              vector of the particle velocities
- rho:              air density [kg/m^3]
- rel_depower:      value between  0.0 and  1.0
- alpha_depower:    depower angle [degrees]
- rel_steering:     value between -1.0 and +1.0

Updates the vector s.forces of the first parameter.
"""
function calc_aero_forces!(s::KPS4_3L, pos, vel)
    r = s.set.radius
    middle_length = s.set.middle_length
    tip_length = s.set.tip_length
    n = s.set.aero_surfaces
    d_s = s.set.minimum_steering_line_distance
    middle_line, left_line, right_line = s.l_tethers
    w = s.set.width

    ρ = s.rho
    # middle_line = 
    # pos_B, pos_C, pos_D: position of the kite particles B, C, and D
    # v_B,   v_C,   v_D:   velocity of the kite particles B, C, and D
    E, C, D = pos[s.num_E], pos[s.num_C], pos[s.num_D]
    v_c, v_d = vel[s.num_C], vel[s.num_D]
     
    # Define functions
    P_c = 0.5 .* (C+D)
    e_y = s.e_y
    e_z = s.e_z
    e_x = s.e_x
    s.e_x .= e_x; s.e_y .= e_y; s.e_z .= e_z # save the kite reference frame in the state

    E_c = E - e_z * s.set.bridle_center_distance # in the aero calculations, E_c is the center of the circle shape on which the kite lies
    F(α) = E_c + e_y*cos(α)*r - e_z*sin(α)*r
    e_r(α) = (E_c - F(α))/norm(E_c-F(α))

    v_cx = dot(v_c, e_x).*e_x
    v_dx = dot(v_d, e_x).*e_x
    v_dy = dot(v_d, e_y).*e_y
    v_dz = dot(v_d, e_z).*e_z
    v_cy = dot(v_c, e_y).*e_y
    v_cz = dot(v_c, e_z).*e_z
    y_lc = norm(C - P_c)
    y_ld = -norm(D - P_c)

    y_l(α) = cos(α) * r
    v_kite(α) = α < π/2 ?
        ((v_cx - v_dx)./(y_lc - y_ld).*(y_l(α) - y_ld) + v_dx) + v_cy + v_cz :
        ((v_cx - v_dx)./(y_lc - y_ld).*(y_l(α) - y_ld) + v_dx) + v_dy + v_dz
    v_a(α) = s.v_wind - v_kite(α)
    v_a_xr(α) = v_a(α) - (dot(v_a(α), cross(e_r(α), e_x))*cross(e_r(α), e_x))

    length(α) = α < π/2 ?
        (tip_length + (middle_length-tip_length)*α*r/(0.5*w)) :
        (tip_length + (middle_length-tip_length)*(π-α)*r/(0.5*w))

    α_l = π/2 - d_s/(2*r) # TODO: move these to outside of function
    α_r = π/2 + d_s/(2*r)

    function d(α)
        if α < α_l
            return middle_line - left_line
        elseif α > α_r
            return middle_line - right_line
        else
            return (-right_line + left_line) / (α_r - α_l) * (α - α_l) + (middle_line - left_line)
        end
    end
    aoa(α) = v_a_xr(α) != [0.0, 0.0, 0.0] ?
        π - acos2(normalize(v_a_xr(α)) ⋅ e_x) + atan(d(α)/length(α)) :
        atan(d(α)/length(α))

    c_l = calc_cl
    c_d = calc_cd

    dL_dα(α) = 0.5*ρ*(norm(v_a_xr(α)))^2*r*length(α)*c_l(aoa(α)) .* normalize(v_a_xr(α) × (e_r(α) × e_x))
    dD_dα(α) = 0.5*ρ*norm(v_a_xr(α))*r*length(α)*c_d(aoa(α)) .* v_a_xr(α) # the sideways drag cannot be calculated with the C_d formula

    # Calculate the integral
    α_0 = pi/2 - w/2/r
    α_middle = pi/2
    dα = (α_middle - α_0) / n
    L_C = sum(dL_dα(α_0 + dα/2 + i*dα) * dα for i in 1:n)
    L_D = sum(dL_dα(pi - (α_0 + dα/2 + i*dα)) * dα for i in 1:n)
    D_C = sum(dD_dα(α_0 + dα/2 + i*dα) * dα for i in 1:n)
    D_D = sum(dD_dα(pi - (α_0 + dα/2 + i*dα)) * dα for i in 1:n)

    s.lift_force .= L_C + L_D
    s.drag_force .= D_C + D_D

    s.forces[s.num_C] .+= (L_C + D_C)
    s.forces[s.num_D] .+= (L_D + D_D)
end

"""
    inner_loop!(s::KPS4_3L, pos, vel, v_wind_gnd, segments, d_tether)

Calculate the forces, acting on all particles.

Output:
- s.forces
- s.v_wind_tether
"""
@inline function inner_loop!(s::KPS4_3L, pos, vel, v_wind_gnd, segments, d_tether)
    for i in eachindex(s.springs)
        p1 = s.springs[i].p1  # First point nr.
        p2 = s.springs[i].p2  # Second point nr.
        height = 0.5 * (pos[p1][3] + pos[p2][3])
        rho = calc_rho(s.am, height)
        @assert height > 0

        s.v_wind_tether .= calc_wind_factor(s.am, height) * v_wind_gnd
        calc_particle_forces!(s, pos[p1], pos[p2], vel[p1], vel[p2], s.springs[i], segments, d_tether, rho, i)
    end
    nothing
end

"""
    loop!(s::KPS4_3L, pos, vel, posd, veld)

Calculate the vectors s.res1 and calculate s.res2 using loops
that iterate over all tether segments. 
"""
function loop!(s::KPS4_3L, pos, vel, posd, veld)
    L_0      = s.l_tethers[1] / s.set.segments
    
    mass_per_meter = s.set.rho_tether * π * (s.set.d_tether/2000.0)^2

    s.res1[1] .= pos[1] # pos = 0
    s.res2[1] .= vel[1] # vel = 0
    particles = s.set.segments + KITE_PARTICLES_3L + 1
    for i in 2:particles
        s.res1[i] .= vel[i] - posd[i] 
    end
    # Compute the masses and forces
    m_tether_particle = mass_per_meter * s.segment_length
    s.masses[s.set.segments+1] = s.set.kcu_mass + 0.5 * m_tether_particle
    # TODO: check if the next two lines are correct
    damping  = s.set.damping / L_0
    c_spring = s.set.c_spring/L_0 
    for i in 1:s.set.segments
        @inbounds s.masses[i] = m_tether_particle
        @inbounds s.springs[i] = SP(s.springs[i].p1, s.springs[i].p2, s.segment_length, c_spring, damping)
    end
    inner_loop!(s, pos, vel, s.v_wind_gnd, s.set.segments, s.set.d_tether/1000.0)
    for i in 2:particles
        s.res2[i] .= veld[i] - (SVector(0, 0, -G_EARTH) - s.forces[i] / s.masses[i])
    end
end

"""
    residual!(res, yd, y::MVector{S, SimFloat}, s::KPS4, time) where S

    N-point tether model, four points for the kite on top:
    Inputs:
        State vector y   = pos1,  pos2, ... , posn, connection_length1-2, vel1,  vel2, . .., veln, connection_vel1-2, length1-3, reel_out_speed1-3
        Derivative   yd  = posd1, posd2, ..., posdn, connection_length1-2, veld1, veld2, ..., veldn, connection_veld1-2, lengthd1-3, reel_out_speedd1-3
        - Without points 1 2 and 3, because they are stationary.
        - With left and right tether points replaced by connection lengths, so they are described by only 1 number instead of 3.
    Output:
    Residual     res = (res1, res2)
        res1 = vel1-posd1,  ..., 
        res2 = veld1-acc1, ..., 
    Will be solved so that res --> 0

    Additional parameters:
    s: Struct with work variables, type KPS4
    S: The dimension of the state vector
The number of the point masses of the model N = S/6, the state of each point 
is represented by two 3 element vectors.
"""
function residual!(res, yd, y::Vector{SimFloat}, s::KPS4_3L, time)
    S = length(y)
    y_ =  MVector{S, SimFloat}(y)
    yd_ =  MVector{S, SimFloat}(yd)
    residual!(res, yd_, y_, s, time)
end
function residual!(res, yd, y::MVector{S, SimFloat}, s::KPS4_3L, time) where S
    T = S-6 # T: three times the number of particles minus 3 origin particles
    num_particles = div(S-6-4, 6) # total number of 3-dimensional particles in y, so excluding 3 stationary points and 2 wire points
    segments = div(T,6) - KITE_PARTICLES_3L + 1 # should equal s.segments
    println("segments\t", segments)
    # Reset the force vector to zero.
    for i in 1:s.num_A
        s.forces[i] .= SVector(0.0, 0, 0)
    end
    # extract the data for the winch simulation
    lengths,  reel_out_speeds  = y[end-5:end-3],  y[end-2:end]
    lengthsd, reel_out_speedsd = yd[end-5:end-3], yd[end-2:end]

    # extract the data of the particles
    y_  = @view y[1:end-6]
    yd_ = @view yd[1:end-6]

    # unpack the vector y
    pos1 = vcat(
        reshape(SVector{3*num_particles}(y_[1:3*num_particles]), Size(3, num_particles)),
        reshape(SVector{2}(y_[num_particles*3+1:num_particles*3+2]),  Size(1, 2)),
    )
    vel1 = vcat(
        reshape(SVector{3*num_particles}(y_[3*num_particles+3:end-2]), Size(3, num_particles)),
        reshape(SVector{2}(y_[end-1:end]),  Size(1, 2)),
    )

    # unpack the vector yd
    posd1 = vcat(
        reshape(SVector{3*num_particles}(yd_[1:3*num_particles]), Size(3, num_particles)),
        reshape(SVector{2}(yd_[num_particles*3+1:num_particles*3+2]),  Size(1, 2)),
    )
    veld1 = vcat(
        reshape(SVector{3*num_particles}(yd_[3*num_particles+3:end-2]),  Size(3, num_particles)),
        reshape(SVector{2}(yd_[end-1:end]),  Size(1, 2)),
    )

    E, C, D = SVector(pos1[end-5]), SVector(pos1[end-4]), SVector(pos1[end-3])
    _, _, e_z = calc_kite_ref_frame!(s, E, C, D)

    # move the right and left tether connection points to the right places
    pos = SVector{s.num_A}(
        if i<=3
            SVector(0.0,0,0)
        elseif 4 <= i <= s.num_E-3 # tether points
            SVector(pos1[i-3])
        elseif i == s.num_E-2 # left tether connection
            C + e_z*pos1[end-1]
        elseif i == s.num_E-1 # right tether connection
            D + e_z*pos1[end]
        elseif i >= s.num_E # kite points
            SVector(pos1[i-5])
        end for i in 1:s.num_A
    )
    vel = SVector{s.num_A}(
        if i<=3
            SVector(0.0,0,0)
        elseif 4 <= i <= s.num_E-3
            SVector(vel1[i-3])
        elseif i == s.num_E-2 # left tether connection, moving on the kite z axis
            C + e_z*vel1[end-1]
        elseif i == s.num_E-1 # right tether connection
            D + e_z*vel1[end]
        elseif i >= s.num_E
            SVector(vel1[i-5])
        end for i in 1:s.num_A
    )
    posd = SVector{s.num_A}(
        if i<=3
            SVector(0.0,0,0)
        elseif 4 <= i <= s.num_E-3
            SVector(posd1[i-3])
        elseif i == s.num_E-2 # left tether connection
            C + e_z*posd1[end-1]
        elseif i == s.num_E-1 # right tether connection
            D + e_z*posd1[end]
        elseif i >= s.num_E
            SVector(posd1[i-5])
        end for i in 1:s.num_A
    )
    veld = SVector{s.num_A}(
        if i<=3
            SVector(0.0,0,0)
        elseif 4 <= i <= s.num_E-3
            SVector(veld1[i-3])
        elseif i == s.num_E-2 # left tether connection
            C + e_z*veld1[end-1]
        elseif i == s.num_E-1 # right tether connection
            D + e_z*veld1[end]
        elseif i >= s.num_E
            SVector(veld1[i-5])
        end for i in 1:s.num_A
    )
    @assert isfinite(pos[4][3])

    # core calculations
    s.l_tethers = lengths
    s.l_connections = connection_lengths
    s.segment_lengths = lengths ./ segments
    calc_aero_forces!(s, pos, vel)
    loop!(s, pos, vel, posd, veld)

    # winch calculations
    res[end-5:end-3] .= lengthsd - reel_out_speeds
    res[end-2:end] .= reel_out_speedsd - calc_acceleration.(s.motors, s.sync_speeds, reel_out_speeds, norm(s.forces[1]), true)

    # copy and flatten result
    for i in 2:div(T,6)+1
        for j in 1:3
            @inbounds res[3*(i-2)+j]              = s.res1[i][j]
            @inbounds res[3*(div(T,6))+3*(i-2)+j] = s.res2[i][j]
        end
    end
    # copy the position vector for easy debugging
    for i in 1:div(T,6)+1
        @inbounds s.pos[i] .= pos[i]
    end
    s.vel_kite .= vel[end-2]
    s.reel_out_speeds = reel_out_speeds

    @assert isfinite(norm(res))
    s.iter += 1

    nothing
end


# =================== getter functions ====================================================

"""
    calc_height(s::KPS4_3L)

Determine the height of the topmost kite particle above ground.
"""
function calc_height(s::KPS4_3L)
    pos_kite(s)[3]
end

"""
    pos_kite(s::KPS4_3L)

Return the position of the kite (top particle).
"""
function pos_kite(s::KPS4_3L)
    s.pos[end-2]
end

"""
    kite_ref_frame(s::KPS4_3L)

Returns a tuple of the x, y, and z vectors of the kite reference frame.
"""
function kite_ref_frame(s::KPS4_3L)
    s.x, s.y, s.z
end

"""
    winch_force(s::KPS4_3L)

Return the absolute value of the force at the winch as calculated during the last timestep. 
"""
function winch_force(s::KPS4_3L) norm(s.last_force) end

# ==================== end of getter functions ================================================

function spring_forces(s::KPS4_3L)
    forces = zeros(SimFloat, s.set.segments+KITE_SPRINGS_3L)
    for i in 1:s.set.segments
        forces[i] =  s.springs[i].c_spring * (norm(s.pos[i+1] - s.pos[i]) - s.segment_length) * s.stiffness_factor
        if forces[i] > 4000.0
            println("Tether raptures for segment $i !")
        end
    end
    for i in 1:KITE_SPRINGS_3L
        p1 = s.springs[i+s.set.segments].p1  # First point nr.
        p2 = s.springs[i+s.set.segments].p2  # Second point nr.
        pos1, pos2 = s.pos[p1], s.pos[p2]
        spring = s.springs[i+s.set.segments]
        l_0 = spring.length # Unstressed length
        k = spring.c_spring * s.stiffness_factor       # Spring constant 
        segment = pos1 - pos2
        norm1 = norm(segment)
        k1 = 0.25 * k # compression stiffness kite segments
        if (norm1 - l_0) > 0.0
            spring_force = k *  (norm1 - l_0) 
        else 
            spring_force = k1 *  (norm1 - l_0)
        end
        forces[i+s.set.segments] = spring_force
        if norm(s.spring_force) > 4000.0
            println("Bridle brakes for spring $i !")
        end
    end
    forces
end

"""
    find_steady_state!(s::KPS4_3L; prn=false, delta = 0.0, stiffness_factor=0.035)

Find an initial equilibrium, based on the inital parameters
`l_tether`, elevation and `reel_out_speeds`.

    X00: parameters that change the shape of the kite system. There are s.set.segments*6+5 params in total

"""
function find_steady_state!(s::KPS4_3L; prn=false, delta = 0.0, stiffness_factor=0.035)
    s.stiffness_factor = stiffness_factor
    res = zeros(MVector{6*(s.num_A-5)+4+6, SimFloat})
    iter = 0

    # helper function for the steady state finder
    function test_initial_condition!(F, x::Vector)
        x1 = copy(x)
        y0, yd0 = init(s, x1)
        residual!(res, yd0, y0, s, 0.0)
        for i in 1:s.set.segments+KITE_PARTICLES_3L-1
            if i != s.set.segments+KITE_PARTICLES_3L-1
                j = i
            else
                j = i + 1
            end
            # copy the x-component of the residual res2 (acceleration)
            F[i]                               = res[1 + 3*(j-1) + 3*(s.set.segments+KITE_PARTICLES_3L)]
            # copy the z-component of the residual res2
            F[i+s.set.segments+KITE_PARTICLES_3L] = res[3 + 3*(j-1) + 3*(s.set.segments+KITE_PARTICLES_3L)]
        end
        # copy the acceleration of point KCU in x direction
        i = s.set.segments+1
        F[end-1]                               = res[1 + 3*(i-1) + 3*(s.set.segments+KITE_PARTICLES_3L)] 
        # copy the acceleration of point C in y direction
        i = s.set.segments+3 
        F[end]                                 = res[2 + 3*(i-1) + 3*(s.set.segments+KITE_PARTICLES_3L)] 
        iter += 1
        return nothing 
    end
    if prn println("\nStarted function test_nlsolve...") end
    X00 = zeros(SimFloat, s.set.segments*6+5)
    results = nlsolve(test_initial_condition!, X00, autoscale=true, xtol=2e-7, ftol=2e-7, iterations=s.set.max_iter)
    if prn println("\nresult: $results") end
    init(s, results.zero)
end
