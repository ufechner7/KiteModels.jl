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
const KITE_PARTICLES = 4
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
    rel_vel::T =           zeros(S, 3)
    half_drag_force::T =   zeros(S, 3)
    forces::SVector{P, KVec3} = zeros(SVector{P, KVec3})
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
    for i in 1:se().segments + KiteModels.KITE_PARTICLES + 1 
        s.forces[i] .= zeros(3)
    end
    s.rho = s.set.rho_0
    s.calc_cl = Spline1D(s.set.alpha_cl, s.set.cl_list)
    s.calc_cd = Spline1D(s.set.alpha_cd, s.set.cd_list) 
end

function KPS4(kcu::KCU)
    s = KPS4{SimFloat, KVec3, kcu.set.segments+KITE_PARTICLES+1, kcu.set.segments+KITE_SPRINGS, SP}()
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

# Calculate the initial conditions y0 and yd0. Tether with the initial elevation angle
# se().elevation, particle zero fixed at origin.
function init(s::KPS4)
    delta = 1e-6
    particles = get_particles(s.set.height_k, s.set.h_bridle, s.set.width, s.set.m_k)
    pos = zeros(SVector{s.set.segments+1+KITE_PARTICLES, KVec3})
    vel = zeros(SVector{s.set.segments+1+KITE_PARTICLES, KVec3})
    sin_el, cos_el = sin(s.set.elevation / 180.0 * pi), cos(s.set.elevation / 180.0 * pi)
    for i in 1:s.set.segments+1
        radius = -(i-1) * (s.set.l_tether/s.set.segments)
        pos[i] .= [-cos_el * radius, delta, -sin_el * radius]
        if i == 1
            vel[i] .= [delta, delta, delta]
        else
            vel[i] .= [delta, delta, 0]
        end
    end
    for i in 1:KITE_PARTICLES
        pos[s.set.segments+1+i] .= particles[i+2]
        vel[s.set.segments+1+i] .= [delta, delta, delta]
    end
    pos, vel
end

""" 
Calculate the drag force of the tether segment, defined by the parameters pos1, pos2, vel1 and vel2
and distribute it equally on the two particles, that are attached to the segment.
The result is stored in the array s.forces. 
"""
function calc_particle_forces(s, pos1, pos2, vel1, vel2, v_wind_tether, spring, stiffnes_factor, segments, d_tether, rho, i)
    l_0 = spring.length # Unstressed length
    k = spring.c_spring * stiffnes_factor       # Spring constant
    c = spring.damping  # Damping coefficient    
    s.segment .= pos1 - pos2
    rel_vel = vel1 - vel2
    s.av_vel .= 0.5 * (vel1 + vel2)
    norm1 = norm(s.segment)
    s.unit_vector .= s.segment / norm1

    k1 = 0.25 * k # compression stiffness kite segments
    k2 = 0.1 * k  # compression stiffness tether segments
    c1 = 6.0 * c  # damping kite segments
    spring_vel   = dot(s.unit_vector, rel_vel)
    if (norm1 - l_0) > 0.0
        if i > segments  # kite springs
             s.spring_force .= (k *  (norm1 - l_0) + (c1 * spring_vel)) * s.unit_vector
        else
             s.spring_force .= (k *  (norm1 - l_0) + (c * spring_vel)) * s.unit_vector
        end
    elseif i > segments # kite spring
        s.spring_force .= (k1 *  (norm1 - l_0) + (c * spring_vel)) * s.unit_vector
    else
        s.spring_force .= (k2 *  (norm1 - l_0) + (c * spring_vel)) * s.unit_vector
    end

    s.v_apparent .= v_wind_tether - s.av_vel
    # TODO: check why d_brindle is not used !!!
    area = norm1 * d_tether
    s.v_app_perp .= s.v_apparent - dot(s.v_apparent, s.unit_vector) * s.unit_vector
    # TODO check the factors 0.25 !!!
    s.half_drag_force .= (-0.25 * rho * s.set.cd_tether * norm(s.v_app_perp) * area) * s.v_app_perp 

    s.forces[spring.p1] .+= s.half_drag_force + s.spring_force
    s.forces[spring.p2] .+= s.half_drag_force - s.spring_force
    nothing
end

function fastlog2(x::Float32)::Float32
    y = Float32(reinterpret(Int32, x))
    y *= 1.1920928955078125f-7
    y - 126.94269504f0
end
function fastlog2(x::Float64)::Float32
   fastlog2(Float32(x))
end

# https://github.com/etheory/fastapprox/blob/master/fastapprox/src/fastexp.h
function fastpow2(x::Float32)::Float32
    clipp = x < -126.0f0 ? -126.0f0 : x
    clipp = min(126f0, max(-126f0, x))
    reinterpret(Float32, UInt32((1 << 23) * (clipp + 126.94269504f0)))
end
function fastpow2(x::Float64)::Float32
   fastpow2(Float32(x))
end

# https://github.com/etheory/fastapprox/blob/master/fastapprox/src/fastpow.h
function fastpow(x::Real, y::Real)::Real
    fastpow2(y * fastlog2(x))
end

"""
Calculate the forces, acting on all particles.
v_wind_tether: out parameter
forces:        out parameter
"""
function inner_loop2(s, pos, vel, v_wind_gnd, stiffnes_factor, segments, d_tether)
    for i in 1:length(s.springs)
        p1 = s.springs[i].p1  # First point nr.
        p2 = s.springs[i].p2  # Second point nr.
        height = 0.5 * (pos[p1][3] + pos[p2][3])
        rho = calc_rho(s, height)
        s.v_wind_tether .= calc_wind_factor(s, height) * v_wind_gnd
        calc_particle_forces(s, pos[p1], pos[p2], vel[p1], vel[p2], s.v_wind_tether, s.springs[i], stiffnes_factor, segments, d_tether, rho, i)
    end
    nothing
end


# 10.8253175473 [ 7.6157172338514298  0.1087959604835919  0.                ] rho:  1.22344998567
# 32.4759526419 [ 8.9098622531401439  0.1272837464734306  0.                ] rho:  1.22035583831
# 54.1265877365 [ 9.584372147186901  0.13691960210267   0.               ] rho:  1.21726951615
# 75.7772228311 [ 10.0563204069071457   0.1436617200986735   0.                ] rho:  1.21419099942
# 97.4278579257 [ 10.4239223873841915   0.1489131769626313   0.                ] rho:  1.21112026835
# 119.07849302 [ 10.7270719657425531   0.1532438852248936   0.                ] rho:  1.20805730327
# 132.273682808 [ 10.8893312909974842   0.1555618755856784   0.                ] rho:  1.20619434991
# 135.361063374 [ 10.9252827098939669   0.1560754672841995   0.                ] rho:  1.20575887521
# 135.112953374 [ 10.9224196754786398   0.156034566792552    0.                ] rho:  1.20579386529
# 135.112953374 [ 10.9224196754786398   0.156034566792552    0.                ] rho:  1.20579386529
# 132.025572808 [ 10.8864110184737015   0.1555201574067672   0.                ] rho:  1.20622935263
# 132.025572808 [ 10.8864110184737015   0.1555201574067672   0.                ] rho:  1.20622935263
# 134.147335048 [ 10.9112339696715832   0.1558747709953083   0.                ] rho:  1.2059300527
# 134.395445048 [ 10.9141146433112723   0.1559159234758753   0.                ] rho:  1.20589505867
# 134.395445048 [ 10.9141146433112723   0.1559159234758753   0.                ] rho:  1.20589505867