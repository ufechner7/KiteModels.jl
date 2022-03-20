# implementation of a four point kite model
# included from KiteModels.jl

# TODO: finish calculation of initial masses

# Array of connections of bridlepoints.
# First point, second point, unstressed length.
const SPRINGS_INPUT = [0.    1.  150.
                       1.    2.   -1.
                       2.    3.   -1.
                       3.    4.   -1.
                       3.    5.   -1.
                       4.    1.   -1.
                       5.    1.   -1.
                       4.    5.   -1.
                       4.    2.   -1.
                       5.    2.   -1.]

# struct, defining the phyical parameters of one spring
@with_kw struct Spring{I, S}
    p1::I = 1         # number of the first point
    p2::I = 2         # number of the second point
    length::S = 1.0   # current unstressed spring length
    c_spring::S = 1.0 # spring constant [N/m]
    damping::S  = 0.1 # damping coefficent [Ns/m]
end

const SP = Spring{Int16, Float64}
const KITE_POINTS = 4
const KITE_SPRINGS = 9
const PRE_STRESS  = 0.9998   # Multiplier for the initial spring lengths.

function zero(::Type{SP})
    SP(0,0,0,0,0)
end

"""
    mutable struct KPS4{S, T, P} <: AbstractKiteModel

State of the kite power system, using a 4 point kite model. Parameters:
- S: Scalar type, e.g. SimFloat
  In the documentation mentioned as Any, but when used in this module it is always SimFloat and not Any.
- T: Vector type, e.g. MVector{3, SimFloat}
- P: number of points of the system, segments+1

Normally a user of this package will not have to access any of the members of this type directly,
use the input and output functions instead.

$(TYPEDFIELDS)
"""
@with_kw mutable struct KPS4{S, T, P, Q, SP} <: AbstractKiteModel
    "Reference to the settings struct"
    set::Settings = se()
    "Reference to the KCU struct (Kite Control Unit, type from the module KitePodSimulor"
    kcu::KCU = KCU()
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
    v_app_perp::T =       zeros(S, 3)
    drag_force::T =       zeros(S, 3)
    lift_force::T =       zeros(S, 3)
    steering_force::T =   zeros(S, 3)
    last_force::T =       zeros(S, 3)
    spring_force::T =     zeros(S, 3)
    total_forces::T =     zeros(S, 3)
    force::T =            zeros(S, 3)
    unit_vector::T =      zeros(S, 3)
    av_vel::T =           zeros(S, 3)
    kite_y::T =           zeros(S, 3)
    segment::T =          zeros(S, 3)
    last_tether_drag::T = zeros(S, 3)
    acc::T =              zeros(S, 3)     
    vec_z::T =            zeros(S, 3)
    pos_kite::T =         zeros(S, 3)
    v_kite::T =           zeros(S, 3)        
    res1::SVector{P, KVec3} = zeros(SVector{P, KVec3})
    res2::SVector{P, KVec3} = zeros(SVector{P, KVec3})
    pos::SVector{P, KVec3} = zeros(SVector{P, KVec3})
    "area of one tether segment"
    seg_area::S =         zero(S) 
    bridle_area::S =      zero(S)
    length::S =           0.0
    area::S =             zero(S)
    last_v_app_norm_tether::S = zero(S)
    "lift coefficient of the kite, depending on the angle of attack"
    param_cl::S =         0.2
    "drag coefficient of the kite, depending on the angle of attack"
    param_cd::S =         1.0
    v_app_norm::S =       zero(S)
    cor_steering::S =     zero(S)
    "azimuth angle in radian; inital value is zero"
    psi::S =              zero(S)
    "elevation angle in radian; initial value about 70 degrees"
    beta::S =             deg2rad(se().elevation)
    last_alpha::S =        0.1
    alpha_depower::S =     0.0
    "relative start time of the current time interval"
    t_0::S =               0.0
    v_reel_out::S =        0.0
    last_v_reel_out::S =   0.0
    l_tether::S =          0.0
    rho::S =               0.0
    depower::S =           0.0
    steering::S =          0.0
    "initial masses of the point masses"
    initial_masses::MVector{P, S} = ones(P)
    "current masses, depending on the total tether length"
    masses::MVector{P, S}         = zeros(P)
    springs::MVector{Q, SP}       = zeros(SP, Q)
end

function clear(s::KPS4)
    s.t_0 = 0.0                              # relative start time of the current time interval
    s.v_reel_out = 0.0
    s.last_v_reel_out = 0.0
    s.area = s.set.area
    s.v_wind        .= [s.set.v_wind, 0, 0]    # wind vector at the height of the kite
    s.v_wind_gnd    .= [s.set.v_wind, 0, 0]    # wind vector at reference height
    s.v_wind_tether .= [s.set.v_wind, 0, 0]
    s.v_apparent    .= [s.set.v_wind, 0, 0]
    s.l_tether = s.set.l_tether
    s.length = s.l_tether / s.set.segments
    s.pos_kite, s.v_kite = zeros(SimFloat, 3), zeros(SimFloat, 3)
    s.beta = deg2rad(s.set.elevation)
    init_masses(s)
    init_springs(s)
    s.rho = s.set.rho_0
    s.calc_cl = Spline1D(s.set.alpha_cl, s.set.cl_list)
    s.calc_cd = Spline1D(s.set.alpha_cd, s.set.cd_list) 
end

function KPS4(kcu::KCU)
    s = KPS4{SimFloat, KVec3, kcu.set.segments+KITE_POINTS+1, kcu.set.segments+KITE_SPRINGS, SP}()
    s.set = kcu.set
    s.kcu = kcu
    s.calc_cl = Spline1D(s.set.alpha_cl, s.set.cl_list)
    s.calc_cd = Spline1D(s.set.alpha_cd, s.set.cd_list)       
    clear(s)
    return s
end

""" 
Calculate the initial orientation of the kite based on the last tether segment and
the apparent wind speed.

Parameters:
vec_c: (pos_n-2) - (pos_n-1) n: number of particles without the three kite particles
                                that do not belong to the main thether (P1, P2 and P3).
Returns:
x, y, z:  the unit vectors of the kite reference frame in the ENU reference frame
"""
function initial_kite_ref_frame(vec_c, v_app)
    z = normalize(vec_c)
    y = normalize(cross(v_app, vec_c))
    x = normalize(cross(y, vec_c))
    return (x, y, z)    
end

""" 
Calculate the initial positions of the particels representing 
a 4-point kite, connected to a kite control unit (KCU). 
"""
function get_particles(height_k, height_b, width, m_k)
    vec_c    = [-15., 0., -25.98076211]
    v_app    = [10.4855, 0, -3.08324]
    pos_pod  = [ 75., 0., 129.90381057]

    # inclination angle of the kite; beta = atan2(-pos_kite[2], pos_kite[1]) ???
    beta = pi/2.0
    x, y, z = initial_kite_ref_frame(vec_c, v_app)

    h_kx = height_k * cos(beta); # print 'h_kx: ', h_kx
    h_kz = height_k * sin(beta); # print 'h_kz: ', h_kz
    h_bx = height_b * cos(beta)
    h_bz = height_b * sin(beta)
    pos_kite = pos_pod - (h_kz + h_bz) * z + (h_kx + h_bx) * x
    pos3 = pos_kite + h_kz * z + 0.5 * width * y + h_kx * x
    pos1 = pos_kite + h_kz * z + (h_kx + width * m_k) * x
    pos4 = pos_kite + h_kz * z - 0.5 * width * y + h_kx * x
    pos0 = pos_kite + (h_kz + h_bz) * z + (h_kx + h_bx) * x
    [zeros(3), pos0, pos1, pos_kite, pos3, pos4]
end

function init_springs(s)
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

function init_masses(s)
    s.masses = zeros(s.set.segments+KITE_POINTS+1)
    l_0 = s.set.l_tether / s.set.segments 
    for i in range(1, s.set.segments)
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

""" 
Calculate the drag force of the tether segment, defined by the parameters pos1, pos2, vel1 and vel2
and distribute it equally on the two particles, that are attached to the segment.
The result is stored in the array s.forces. 
"""
function calc_particle_forces(s, pos1, pos2, vel1, vel2, v_wind_tether, spring, forces, stiffnes_factor, segments, d_tether, i)
    p_1 = spring.p1     # Index of point nr. 1
    p_2 = spring.p2     # Index of point nr. 2
    l_0 = spring.length # Unstressed length
    k = spring.c_spring * stiffnes_factor       # Spring constant
    c = spring.damping  # Damping coefficient    
    s.segment .= pos1 - pos2
    rel_vel = vel1 - vel2
    s.av_vel .= 0.5 * (vel1 + vel2)
    norm1 = norm(s.segment)
    unit_vector = s.segment / norm1
    k1 = 0.25 * k # compression stiffness kite segments
    k2 = 0.1 * k  # compression stiffness tether segments
    spring_vel   = dot(unit_vector, rel_vel)
end

# def calcParticleForces_(pos1, pos2, vel1, vel2, v_wind_tether, spring, forces, stiffnes_factor, segments, \
#                        d_tether, i):
#     p_1 = int_(spring[0])     # Index of point nr. 1
#     p_2 = int_(spring[1])     # Index of point nr. 2
#     l_0 = spring[2]     # Unstressed length
#     k = spring[3] * stiffnes_factor       # Spring constant
#     c = spring[4]       # Damping coefficient
#     segment = pos1 - pos2
#     rel_vel = vel1 - vel2
#     av_vel = 0.5 * (vel1 + vel2)
#     norm1 = la.norm(segment)
#     unit_vector = segment / norm1
#     k1 = 0.25 * k # compression stiffness kite segments
#     k2 = 0.1 * k # compression stiffness tether segments
#     spring_vel   = la.dot(unit_vector, rel_vel)
#     if (norm1 - l_0) > 0.0:
#         spring_force = (k *  (norm1 - l_0) + (c * spring_vel)) * unit_vector
#     elif i >= segments: # kite spring
#         spring_force = (k1 *  (norm1 - l_0) + (c * spring_vel)) * unit_vector
#     else:
#         spring_force = (k2 *  (norm1 - l_0) + (c * spring_vel)) * unit_vector
#     # Aerodynamic damping for particles of the tether and kite
#     v_apparent = v_wind_tether - av_vel
#     area = norm1 * d_tether
#     v_app_perp = v_apparent - la.dot(v_apparent, unit_vector) *unit_vector
#     half_drag_force = -0.25 * 1.25 * C_D_TETHER * la.norm(v_app_perp) * area * v_app_perp

#     forces[p_1] += half_drag_force + spring_force
#     forces[p_2] += half_drag_force - spring_force
