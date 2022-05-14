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
configured in the file data/settings.yaml). The kite is modelled using 4 point masses and 3 aerodynamic 
surfaces. The spring constant and the damping decrease with the segment length. The aerodynamic kite forces
are acting on three of the four kite point masses. 

Four point kite model, included from KiteModels.jl.

Scientific background: http://arxiv.org/abs/1406.6218 =#

# Array of connections of bridlepoints.
# First point, second point, unstressed length.
const SPRINGS_INPUT = [0.    1.  150.
                       1.    2.   -1. # s1, p7, p8
                       4.    2.   -1. # s2, p10, p8                        
                       4.    5.   -1. # s3, p10, p11
                       3.    4.   -1. # s4, p9, p10
                       5.    1.   -1. # s5, p11, p7
                       4.    1.   -1. # s6, p10, p7
                       3.    5.   -1. # s7, p9, p11
                       5.    2.   -1. # s8, p11, p8
                       2.    3.   -1.] # s9, p8, p9

# struct, defining the phyical parameters of one spring
@with_kw struct Spring{I, S}
    p1::I = 1         # number of the first point
    p2::I = 2         # number of the second point
    length::S = 1.0   # current unstressed spring length
    c_spring::S = 1.0 # spring constant [N/m]
    damping::S  = 0.1 # damping coefficent [Ns/m]
end

const SP = Spring{Int16, Float64}
const KITE_PARTICLES = 4
const KITE_SPRINGS = 9
const KITE_ANGLE = 3.83 # angle between the kite and the last tether segment due to the mass of the control pod
const PRE_STRESS  = 0.9998   # Multiplier for the initial spring lengths.
const KS = deg2rad(16.565 * 1.064 * 0.875 * 1.033 * 0.9757 * 1.083)  # max steering
const DRAG_CORR = 0.93       # correction of the drag for the 4-point model
function zero(::Type{SP})
    SP(0,0,0,0,0)
end

"""
    mutable struct KPS4{S, T, P, Q, SP} <: AbstractKiteModel

State of the kite power system, using a 4 point kite model. Parameters:
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
@with_kw mutable struct KPS4{S, T, P, Q, SP} <: AbstractKiteModel
    "Reference to the settings struct"
    set::Settings = se()
    "Reference to the KCU model (Kite Control Unit as implemented in the package KitePodModels"
    kcu::KCU = KCU()
    "Reference to the atmospheric model as implemented in the package AtmosphericModels"
    am::AtmosphericModel = AtmosphericModel()
    "Reference to winch model as implemented in the package WinchModels"
    wm::AbstractWinchModel = AsyncMachine()
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
    segment::T =          zeros(S, 3) 
    "a copy of the residual one (pos,vel) for debugging and unit tests"    
    res1::SVector{P, KVec3} = zeros(SVector{P, KVec3})
    "a copy of the residual two (vel,acc) for debugging and unit tests"
    res2::SVector{P, KVec3} = zeros(SVector{P, KVec3})
    "a copy of the actual positions as output for the user"
    pos::SVector{P, KVec3} = zeros(SVector{P, KVec3})
    "velocity vector of the kite"
    vel_kite::T =          zeros(S, 3)
    "unstressed segment length [m]"
    segment_length::S =           0.0
    "lift coefficient of the kite, depending on the angle of attack"
    param_cl::S =         0.2
    "drag coefficient of the kite, depending on the angle of attack"
    param_cd::S =         1.0
    "azimuth angle in radian; inital value is zero"
    psi::S =              zero(S)
    "elevation angle in radian; initial value about 70 degrees"
    beta::S =             deg2rad(se().elevation)
    alpha_depower::S =     0.0
    "relative start time of the current time interval"
    t_0::S =               0.0
    "reel out speed of the winch"
    v_reel_out::S =        0.0
    "reel out speed at the last time step"
    last_v_reel_out::S =   0.0
    "unstretched tether length"
    l_tether::S =          0.0
    "air density"
    rho::S =               0.0
    "actual relative depower setting,  must be between    0 .. 1.0"
    depower::S =           0.0
    "actual relative steering setting, must be between -1.0 .. 1.0"
    steering::S =          0.0
    "multiplier for the stiffniss of tether and bridle"
    stiffness_factor::S =  1.0
    "initial masses of the point masses"
    initial_masses::MVector{P, S} = ones(P)
    "current masses, depending on the total tether length"
    masses::MVector{P, S}         = zeros(P)
    springs::MVector{Q, SP}       = zeros(SP, Q)
    rel_vel::T =           zeros(S, 3)
    half_drag_force::T =   zeros(S, 3)
    forces::SVector{P, KVec3} = zeros(SVector{P, KVec3})
    "synchronous speed of the motor/ generator"
    sync_speed::S =        0.0         
    x::T =                 zeros(S, 3)
    y::T =                 zeros(S, 3)
    z::T =                 zeros(S, 3)
end

"""
    clear!(s::KPS4)

Initialize the kite power model.
"""
function clear!(s::KPS4)
    s.t_0 = 0.0                              # relative start time of the current time interval
    s.v_reel_out = 0.0
    s.last_v_reel_out = 0.0
    s.v_wind_gnd    .= [s.set.v_wind, 0, 0]    # wind vector at reference height
    s.v_wind_tether .= [s.set.v_wind, 0, 0]
    s.v_apparent    .= [s.set.v_wind, 0, 0]
    height = sin(deg2rad(s.set.elevation)) * (s.set.l_tether)
    s.v_wind .= s.v_wind_gnd * calc_wind_factor(s.am, height)

    s.l_tether = s.set.l_tether
    s.segment_length = s.l_tether / s.set.segments
    s.beta = deg2rad(s.set.elevation)
    init_masses!(s)
    init_springs!(s)
    for i in 1:se().segments + KiteModels.KITE_PARTICLES + 1 
        s.forces[i] .= zeros(3)
    end
    s.drag_force .= [0.0, 0, 0]
    s.lift_force .= [0.0, 0, 0]
    s.rho = s.set.rho_0
    s.calc_cl = Spline1D(s.set.alpha_cl, s.set.cl_list)
    s.calc_cd = Spline1D(s.set.alpha_cd, s.set.cd_list) 
    s.kcu.depower = s.set.depower/100.0
    s.kcu.set_depower = s.kcu.depower
    KiteModels.set_depower_steering!(s, get_depower(s.kcu), get_steering(s.kcu))
end

function KPS4(kcu::KCU)
    s = KPS4{SimFloat, KVec3, kcu.set.segments+KITE_PARTICLES+1, kcu.set.segments+KITE_SPRINGS, SP}()
    s.set = kcu.set
    s.kcu = kcu
    s.calc_cl = Spline1D(s.set.alpha_cl, s.set.cl_list)
    s.calc_cd = Spline1D(s.set.alpha_cd, s.set.cd_list)       
    clear!(s)
    return s
end

"""
    initial_kite_ref_frame(vec_c, v_app)

Calculate the initial orientation of the kite based on the last tether segment and
the apparent wind speed.

Parameters:
- vec_c: (pos_n-2) - (pos_n-1) n: number of particles without the three kite particles
                                  that do not belong to the main thether (P1, P2 and P3).
- v_app: vector of the apparent wind speed

Returns:
x, y, z:  the unit vectors of the kite reference frame in the ENU reference frame
"""
function initial_kite_ref_frame(vec_c, v_app)
    z = normalize(vec_c)
    y = normalize(v_app × vec_c)
    x = normalize(y × vec_c)
    x, y, z
end

"""
    get_particles(height_k, height_b, width, m_k, pos_pod= [ 75., 0., 129.90381057], vec_c=[-15., 0., -25.98076211], v_app=[10.4855, 0, -3.08324])

Calculate the initial positions of the particels representing 
a 4-point kite, connected to a kite control unit (KCU). 

Parameters:
- height_k: height of the kite itself, not above ground [m]
- height_b: height of the bridle [m]
- width: width of the kite [m]
- mk: relative nose distance
- pos_pod: position of the control pod
- vec_c: vector of the last tether segment
"""
function get_particles(height_k, height_b, width, m_k, pos_pod= [ 75., 0., 129.90381057], vec_c=[-15., 0., -25.98076211], v_app=[10.4855, 0, -3.08324])
    # inclination angle of the kite; beta = atan(-pos_kite[2], pos_kite[1]) ???
    beta = pi/2.0
    x, y, z = initial_kite_ref_frame(vec_c, v_app)

    h_kx = height_k * cos(beta); # print 'h_kx: ', h_kx
    h_kz = height_k * sin(beta); # print 'h_kz: ', h_kz
    h_bx = height_b * cos(beta)
    h_bz = height_b * sin(beta)
    pos_kite = pos_pod - (h_kz + h_bz) * z + (h_kx + h_bx) * x   # top,        poing B in diagram
    pos_C = pos_kite + h_kz * z + 0.5 * width * y + h_kx * x     # side point, point C in diagram
    pos_A = pos_kite + h_kz * z + (h_kx + width * m_k) * x       # nose,       point A in diagram
    pos_D = pos_kite + h_kz * z - 0.5 * width * y + h_kx * x     # side point, point D in diagram
    pos0 = pos_kite + (h_kz + h_bz) * z + (h_kx + h_bx) * x      # equal to pos_pod, P_KCU in diagram
    [zeros(3), pos0, pos_A, pos_kite, pos_C, pos_D] # 0, p7, p8, p9, p10, p11
end

"""
    calc_height(s::KPS4)

Determine the height of the topmost kite particle above ground.
"""
function calc_height(s::KPS4)
    pos_kite(s)[3]
end

"""
    pos_kite(s::KPS4)

Return the position of the kite (top particle).
"""
function pos_kite(s::KPS4)
    s.pos[end-2]
end

function init_springs!(s::KPS4)
    l_0     = s.set.l_tether / s.set.segments 
    particles = get_particles(s.set.height_k, s.set.h_bridle, s.set.width, s.set.m_k)
    for j in 1:size(SPRINGS_INPUT)[1]
        # build the tether segments
        if j == 1
            for i in 1:s.set.segments
                k = s.set.e_tether * (s.set.d_tether/2000.0)^2 * pi  / l_0  # Spring stiffness for this spring [N/m]
                c = s.set.damping/l_0                                       # Damping coefficient [Ns/m]
                s.springs[i] = SP(i, i+1, l_0, k, c)
            end
        # build the bridle segments
        else
            p0, p1 = SPRINGS_INPUT[j, 1]+1, SPRINGS_INPUT[j, 2]+1 # point 0 and 1
            if SPRINGS_INPUT[j, 3] == -1
                l_0 = norm(particles[Int(p1)] - particles[Int(p0)]) * PRE_STRESS
                k = s.set.e_tether * (s.set.d_line/2000.0)^2 * pi / l_0
                p0 += s.set.segments - 1 # correct the index for the start and end particles of the bridle
                p1 += s.set.segments - 1
                c = s.set.damping/ l_0
                s.springs[j+s.set.segments-1] = SP(Int(p0), Int(p1), l_0, k, c)
            end
        end
    end
    s.springs
end

function init_masses!(s::KPS4)
    s.masses = zeros(s.set.segments+KITE_PARTICLES+1)
    l_0 = s.set.l_tether / s.set.segments 
    for i in 1:s.set.segments
        s.masses[i]   += 0.5 * l_0 * s.set.rho_tether * (s.set.d_tether/2000.0)^2 * pi
        s.masses[i+1] += 0.5 * l_0 * s.set.rho_tether * (s.set.d_tether/2000.0)^2 * pi
    end
    s.masses[s.set.segments+1] += s.set.kcu_mass
    k2 = s.set.rel_top_mass * (1.0 - s.set.rel_nose_mass)
    k3 = 0.5 * (1.0 - s.set.rel_top_mass) * (1.0 - s.set.rel_nose_mass)
    k4 = 0.5 * (1.0 - s.set.rel_top_mass) * (1.0 - s.set.rel_nose_mass)
    s.masses[s.set.segments+2] += s.set.rel_nose_mass * s.set.mass
    s.masses[s.set.segments+3] += k2 * s.set.mass
    s.masses[s.set.segments+4] += k3 * s.set.mass
    s.masses[s.set.segments+5] += k4 * s.set.mass  
    s.masses 
end

# Calculate the initial vectors pos and vel. Tether with the initial elevation angle
# se().elevation, particle zero fixed at origin.
# X is a vector of deviations in x and z positions, to be varied to find the inital equilibrium
function init_pos_vel(s::KPS4, X=zeros(2 * (s.set.segments+KITE_PARTICLES)))
    pos, vel, acc = init_pos_vel_acc(s, X; old=true)
    pos, vel
end

function init_pos_vel_acc(s::KPS4, X=zeros(2 * (s.set.segments+KITE_PARTICLES)+1); old=false, delta = 0.0)
    pos = zeros(SVector{s.set.segments+1+KITE_PARTICLES, KVec3})
    vel = zeros(SVector{s.set.segments+1+KITE_PARTICLES, KVec3})
    acc = zeros(SVector{s.set.segments+1+KITE_PARTICLES, KVec3})
    pos[1] .= [0.0, delta, 0.0]
    vel[1] .= [delta, delta, delta]
    acc[1] .= [delta, delta, delta]
    sin_el, cos_el = sin(s.set.elevation / 180.0 * pi), cos(s.set.elevation / 180.0 * pi)
    for i in 1:s.set.segments
        radius = -i * (s.set.l_tether/s.set.segments)
        pos[i+1] .= [-cos_el * radius + X[i], delta, -sin_el * radius + X[s.set.segments+KITE_PARTICLES-1+i]]
        vel[i+1] .= [delta, delta, 0]
        acc[i+1] .= [delta, delta, -9.81]
    end
    vec_c = pos[s.set.segments] - pos[s.set.segments+1]
    if old
        particles = get_particles(s.set.height_k, s.set.h_bridle, s.set.width, s.set.m_k)
    else
        particles = get_particles(s.set.height_k, s.set.h_bridle, s.set.width, s.set.m_k, pos[s.set.segments+1], rotate_in_xz(vec_c, deg2rad(KITE_ANGLE)), s.v_apparent)
    end
    j = 1
    for i in [1,2,3] # set p8, p9, p10
        pos[s.set.segments+1+i] .= particles[i+2] + [X[s.set.segments+j], 0, X[2*s.set.segments+KITE_PARTICLES-1+j]]
        vel[s.set.segments+1+i] .= [delta, delta, delta]
        acc[s.set.segments+1+i] .= [delta, delta, -9.81]
        j +=1
    end
    acc[s.set.segments+1+4] .= [delta, delta, -9.81]
    vel[s.set.segments+1+4] .= [delta, delta, delta]
    # set p10=C and p11=D
    # x and z component of the right and left particle must be equal
    pos[s.set.segments+1+4][1] = pos[s.set.segments+1+3][1]  # D.x = C.x
    pos[s.set.segments+1+4][3] = pos[s.set.segments+1+3][3]  # D.z = C.z
    pos[s.set.segments+1+3][2] += X[end]                     # Y position of point C
    pos[s.set.segments+1+4][2] = -pos[s.set.segments+1+3][2] # Y position of point D
    for i in 1:length(pos)
        s.pos[i] .= pos[i]
    end
    pos, vel, acc
end

function init_inner(s::KPS4, X=zeros(2 * (s.set.segments+KITE_PARTICLES-1)+1); old=false, delta=0.0)
    pos, vel, acc = init_pos_vel_acc(s, X; old=old, delta=delta)
    vcat(pos[2:end], vel[2:end]), vcat(vel[2:end], acc[2:end])
end

# same as above, but returns a tuple of two one dimensional arrays
function init(s::KPS4, X=zeros(2 * (s.set.segments+KITE_PARTICLES-1)+1); old=false, delta=0.0)
    res1_, res2_ = init_inner(s, X; old=old, delta = delta)
    res1, res2  = vcat(reduce(vcat, res1_), [s.l_tether, 0]), vcat(reduce(vcat, res2_),[0,0])
    MVector{6*(s.set.segments+KITE_PARTICLES)+2, Float64}(res1), MVector{6*(s.set.segments+KITE_PARTICLES)+2, Float64}(res2)
end

""" 
    calc_particle_forces!(s::KPS4, pos1, pos2, vel1, vel2, spring, segments, d_tether, rho, i)

Calculate the drag force of the tether segment, defined by the parameters pos1, pos2, vel1 and vel2
and distribute it equally on the two particles, that are attached to the segment.
The result is stored in the array s.forces. 
"""
@inline function calc_particle_forces!(s, pos1, pos2, vel1, vel2, spring, segments, d_tether, rho, i)
    l_0 = spring.length # Unstressed length
    k = spring.c_spring * s.stiffness_factor  # Spring constant
    c = spring.damping                        # Damping coefficient    
    s.segment .= pos1 - pos2
    rel_vel = vel1 - vel2
    av_vel = 0.5 * (vel1 + vel2)
    norm1 = norm(s.segment)
    unit_vector = s.segment / norm1

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
    # TODO: check why d_brindle is not used !!!
    area = norm1 * d_tether
    v_app_perp = s.v_apparent - s.v_apparent ⋅ unit_vector * unit_vector
    s.half_drag_force .= (-0.25 * rho * s.set.cd_tether * norm(v_app_perp) * area) * v_app_perp 

    @inbounds s.forces[spring.p1] .+= s.half_drag_force + s.spring_force
    @inbounds s.forces[spring.p2] .+= s.half_drag_force - s.spring_force
    if i == 1 s.last_force .= s.forces[spring.p1] end
    nothing
end

"""
    inner_loop!(s::KPS4, pos, vel, v_wind_gnd, segments, d_tether)

Calculate the forces, acting on all particles.

Output:
- s.forces
- `s.v_wind_tether`
"""
@inline function inner_loop!(s::KPS4, pos, vel, v_wind_gnd, segments, d_tether)
    for i in 1:length(s.springs)
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
    winch_force(s::KPS4)

Return the absolute value of the force at the winch as calculated during the last timestep. 
"""
function winch_force(s::KPS4) norm(s.last_force) end

"""
    calc_aero_forces!(s::KPS4, pos, vel, rho, alpha_depower, rel_steering)

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
function calc_aero_forces!(s::KPS4, pos, vel, rho, alpha_depower, rel_steering)
    rel_side_area = s.set.rel_side_area/100.0    # defined in percent
    K = 1 - rel_side_area                        # correction factor for the drag
    # pos_B, pos_C, pos_D: position of the kite particles B, C, and D
    # v_B,   v_C,   v_D:   velocity of the kite particles B, C, and D
    pos_B, pos_C, pos_D = pos[s.set.segments+3], pos[s.set.segments+4], pos[s.set.segments+5]
    v_B,   v_C,   v_D   = vel[s.set.segments+3], vel[s.set.segments+4], vel[s.set.segments+5]
    va_2,  va_3,  va_4  = s.v_wind - v_B, s.v_wind - v_C, s.v_wind - v_D
 
    pos_centre = 0.5 * (pos_C + pos_D)
    delta = pos_B - pos_centre
    z = -normalize(delta)
    y = normalize(pos_C - pos_D)
    x = y × z
    s.x .= x; s.y .= y; s.z .= z # save the kite reference frame in the state

    va_xz2 = va_2 - (va_2 ⋅ y) * y
    va_xy3 = va_3 - (va_3 ⋅ z) * z
    va_xy4 = va_4 - (va_4 ⋅ z) * z

    alpha_2 = rad2deg(π - acos2(normalize(va_xz2) ⋅ x) - alpha_depower)     + s.set.alpha_zero
    alpha_3 = rad2deg(π - acos2(normalize(va_xy3) ⋅ x) - rel_steering * KS) + s.set.alpha_ztip
    alpha_4 = rad2deg(π - acos2(normalize(va_xy4) ⋅ x) + rel_steering * KS) + s.set.alpha_ztip

    CL2, CD2 = calc_cl(alpha_2), DRAG_CORR * calc_cd(alpha_2)
    CL3, CD3 = calc_cl(alpha_3), DRAG_CORR * calc_cd(alpha_3)
    CL4, CD4 = calc_cl(alpha_4), DRAG_CORR * calc_cd(alpha_4)
    L2 = (-0.5 * rho * (norm(va_xz2))^2 * s.set.area * CL2) * normalize(va_2 × y)
    L3 = (-0.5 * rho * (norm(va_xy3))^2 * s.set.area * rel_side_area * CL3) * normalize(va_3 × z)
    L4 = (-0.5 * rho * (norm(va_xy4))^2 * s.set.area * rel_side_area * CL4) * normalize(z × va_4)
    D2 = (-0.5 * K * rho * norm(va_2) * s.set.area * CD2) * va_2
    D3 = (-0.5 * K * rho * norm(va_3) * s.set.area * rel_side_area * CD3) * va_3
    D4 = (-0.5 * K * rho * norm(va_4) * s.set.area * rel_side_area * CD4) * va_4
    s.lift_force .= L2
    s.drag_force .= D2 + D3 + D4

    s.forces[s.set.segments + 3] .+= (L2 + D2)
    s.forces[s.set.segments + 4] .+= (L3 + D3)
    s.forces[s.set.segments + 5] .+= (L4 + D4)
end

"""
    kite_ref_frame(s::KPS4)

Returns a tuple of the x, y, and z vectors of the kite reference frame.
"""
function kite_ref_frame(s::KPS4)
    s.x, s.y, s.z
end

"""
    loop!(s::KPS4, pos, vel, posd, veld)

Calculate the vectors s.res1 and calculate s.res2 using loops
that iterate over all tether segments. 
"""
function loop!(s::KPS4, pos, vel, posd, veld)
    L_0      = s.l_tether / s.set.segments
    # mass_per_meter = s.set.rho_tether * π * (s.set.d_tether/2000.0)^2
    mass_per_meter = 0.011
    s.res1[1] .= pos[1]
    s.res2[1] .= vel[1]
    particles = s.set.segments + KITE_PARTICLES + 1
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
    State vector y   = pos1,  pos2, ... , posn,  vel1,  vel2, . .., veln,  length, v_reel_out
    Derivative   yd  = posd1, posd2, ..., posdn, veld1, veld2, ..., veldn, lengthd, v_reel_outd
    Output:
    Residual     res = res1, res2 = vel1-posd1,  ..., veld1-acc1, ..., 

    Additional parameters:
    s: Struct with work variables, type KPS4
    S: The dimension of the state vector
The number of the point masses of the model N = S/6, the state of each point 
is represented by two 3 element vectors.
"""
function residual!(res, yd, y::MVector{S, SimFloat}, s::KPS4, time) where S
    T = S-2 # T: three times the number of particles excluding the origin
    segments = div(T,6) - KITE_PARTICLES
    
    # Reset the force vector to zero.
    for i in 1:segments+KITE_PARTICLES+1
        s.forces[i] .= SVector(0.0, 0, 0)
    end
    # extract the data for the winch simulation
    length,  v_reel_out  = y[end-1],  y[end]
    lengthd, v_reel_outd = yd[end-1], yd[end]
    # extract the data of the particles
    y_  = @view y[1:end-2]
    yd_ = @view yd[1:end-2]
    # unpack the vectors y and yd
    part  = reshape(SVector{T}(y_),  Size(3, div(T,6), 2))
    partd = reshape(SVector{T}(yd_), Size(3, div(T,6), 2))
    pos1, vel1 = part[:,:,1], part[:,:,2]
    pos = SVector{div(T,6)+1}(if i==1 SVector(0.0,0,0) else SVector(pos1[:,i-1]) end for i in 1:div(T,6)+1)
    vel = SVector{div(T,6)+1}(if i==1 SVector(0.0,0,0) else SVector(vel1[:,i-1]) end for i in 1:div(T,6)+1)
    posd1, veld1 = partd[:,:,1], partd[:,:,2]
    posd = SVector{div(T,6)+1}(if i==1 SVector(0.0,0,0) else SVector(posd1[:,i-1]) end for i in 1:div(T,6)+1)
    veld = SVector{div(T,6)+1}(if i==1 SVector(0.0,0,0) else SVector(veld1[:,i-1]) end for i in 1:div(T,6)+1)
    @assert ! isnan(pos[2][3])

    # core calculations
    s.l_tether = length
    s.segment_length = length / segments
    calc_aero_forces!(s, pos, vel, s.rho, s.alpha_depower, s.steering)
    loop!(s, pos, vel, posd, veld)

    # winch calculations
    res[end-1] = lengthd - v_reel_out
    res[end] = v_reel_outd - calc_acceleration(s.wm, s.sync_speed, v_reel_out, norm(s.forces[1]), true)

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

    @assert ! isnan(norm(res))
    s.iter += 1

    nothing
end

function spring_forces(s::KPS4)
    forces = zeros(SimFloat, s.set.segments+KITE_SPRINGS)
    for i in 1:s.set.segments
        forces[i] =  s.springs[i].c_spring * (norm(s.pos[i+1] - s.pos[i]) - s.segment_length) * s.stiffness_factor
        if forces[i] > 4000.0
            println("Tether raptures for segment $i !")
        end
    end
    for i in 1:KITE_SPRINGS
        p1 = s.springs[i+s.set.segments].p1  # First point nr.
        p2 = s.springs[i+s.set.segments].p2  # Second point nr.
        pos1, pos2 = s.pos[p1], s.pos[p2]
        spring = s.springs[i+s.set.segments]
        l_0 = spring.length # Unstressed length
        k = spring.c_spring * s.stiffness_factor       # Spring constant 
        s.segment .= pos1 - pos2
        norm1 = norm(s.segment)
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
    find_steady_state!(s::KPS4; prn=false, delta = 0.0, stiffness_factor=0.035)

Find an initial equilibrium, based on the inital parameters
`l_tether`, elevation and `v_reel_out`.
"""
function find_steady_state!(s::KPS4; prn=false, delta = 0.0, stiffness_factor=0.035)
    s.stiffness_factor = stiffness_factor
    res = zeros(MVector{6*(s.set.segments+KITE_PARTICLES)+2, SimFloat})
    iter = 0

    # helper function for the steady state finder
    function test_initial_condition!(F, x::Vector)
        x1 = copy(x)
        y0, yd0 = init(s, x1)
        residual!(res, yd0, y0, s, 0.0)
        for i in 1:s.set.segments+KITE_PARTICLES-1
            if i != s.set.segments+KITE_PARTICLES-1
                j = i
            else
                j = i + 1
            end
            # copy the x-component of the residual res2 (acceleration)
            F[i]                               = res[1 + 3*(j-1) + 3*(s.set.segments+KITE_PARTICLES)]
            # copy the z-component of the residual res2
            F[i+s.set.segments+KITE_PARTICLES] = res[3 + 3*(j-1) + 3*(s.set.segments+KITE_PARTICLES)]
        end
        # copy the acceleration of point KCU in x direction
        i = s.set.segments+1
        F[end-1]                               = res[1 + 3*(i-1) + 3*(s.set.segments+KITE_PARTICLES)] 
        # copy the acceleration of point C in y direction
        i = s.set.segments+3 
        F[end]                                 = res[2 + 3*(i-1) + 3*(s.set.segments+KITE_PARTICLES)] 
        iter += 1
        return nothing 
    end
    if prn println("\nStarted function test_nlsolve...") end
    X00 = zeros(SimFloat, 2*(s.set.segments+KITE_PARTICLES-1)+2)
    results = nlsolve(test_initial_condition!, X00, autoscale=true, xtol=2e-7, ftol=2e-7, iterations=s.set.max_iter)
    if prn println("\nresult: $results") end
    init(s, results.zero)
end

# rotate a 3d vector around the y axis
function rotate_in_xz(vec, angle)
    result = similar(vec)
    result[1] = cos(angle) * vec[1] - sin(angle) * vec[3]
    result[2] = vec[2]
    result[3] = cos(angle) * vec[3] + sin(angle) * vec[1]
    result
end
