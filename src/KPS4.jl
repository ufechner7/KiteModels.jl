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

This model implements a 3D mass-spring system with reel-out. It uses five tether segments (the number can be
configured in the file data/settings.yaml). The kite is modelled as additional mass at the end of the tether.
The spring constant and the damping decrease with the segment length. The aerodynamic kite forces are
calculated, depending on reel-out speed, depower and steering settings. 

Scientific background: http://arxiv.org/abs/1406.6218 =#

# implementation of a four point kite model
# included from KiteModels.jl

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
    "spring force of the current tether segment, output of calc_particle_forces"
    spring_force::T =     zeros(S, 3)
    segment::T =          zeros(S, 3)
    v_kite::T =           zeros(S, 3)        
    res1::SVector{P, KVec3} = zeros(SVector{P, KVec3})
    res2::SVector{P, KVec3} = zeros(SVector{P, KVec3})
    pos::SVector{P, KVec3} = zeros(SVector{P, KVec3})
    "unstressed reelout length [m]"
    length::S =           0.0
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
    s.v_wind        .= [s.set.v_wind, 0, 0]    # wind vector at the height of the kite
    s.v_wind_gnd    .= [s.set.v_wind, 0, 0]    # wind vector at reference height
    s.v_wind_tether .= [s.set.v_wind, 0, 0]
    s.v_apparent    .= [s.set.v_wind, 0, 0]
    s.l_tether = s.set.l_tether
    s.length = s.l_tether / s.set.segments
    s.v_kite = zeros(SimFloat, 3)
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
# X is a vector of deviations in x and z positions, to be varied to find the inital equilibrium
function init(s::KPS4, X=zeros(2 * (s.set.segments+KITE_PARTICLES)))
    delta = 1e-6
    particles = get_particles(s.set.height_k, s.set.h_bridle, s.set.width, s.set.m_k)
    pos = zeros(SVector{s.set.segments+1+KITE_PARTICLES, KVec3})
    vel = zeros(SVector{s.set.segments+1+KITE_PARTICLES, KVec3})
    pos[1] .= [0.0, delta, 0.0]
    vel[1] .= [delta, delta, delta]
    sin_el, cos_el = sin(s.set.elevation / 180.0 * pi), cos(s.set.elevation / 180.0 * pi)
    for i in 1:s.set.segments
        radius = -i * (s.set.l_tether/s.set.segments)
        pos[i+1] .= [-cos_el * radius + X[i], delta, -sin_el * radius + X[s.set.segments+i]]
        vel[i+1] .= [delta, delta, 0]
    end
    for i in 1:KITE_PARTICLES
        pos[s.set.segments+1+i] .= particles[i+2] + [X[2*s.set.segments+i], 0, X[2*s.set.segments+KITE_PARTICLES+i]]
        vel[s.set.segments+1+i] .= [delta, delta, delta]
    end
    pos, vel
end

""" 
Calculate the drag force of the tether segment, defined by the parameters pos1, pos2, vel1 and vel2
and distribute it equally on the two particles, that are attached to the segment.
The result is stored in the array s.forces. 
"""
@inline function calc_particle_forces(s, pos1, pos2, vel1, vel2, spring, stiffnes_factor, segments, d_tether, rho, i)
    l_0 = spring.length # Unstressed length
    k = spring.c_spring * stiffnes_factor       # Spring constant
    c = spring.damping  # Damping coefficient    
    s.segment .= pos1 - pos2
    rel_vel = vel1 - vel2
    av_vel = 0.5 * (vel1 + vel2)
    norm1 = norm(s.segment)
    unit_vector = s.segment / norm1

    k1 = 0.25 * k # compression stiffness kite segments
    k2 = 0.1 * k  # compression stiffness tether segments
    c1 = 6.0 * c  # damping kite segments
    spring_vel   = dot(unit_vector, rel_vel)
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
    v_app_perp = s.v_apparent - dot(s.v_apparent, unit_vector) * unit_vector
    # TODO check the factors 0.25 !!!
    s.half_drag_force .= (-0.25 * rho * s.set.cd_tether * norm(v_app_perp) * area) * v_app_perp 

    @inbounds s.forces[spring.p1] .+= s.half_drag_force + s.spring_force
    @inbounds s.forces[spring.p2] .+= s.half_drag_force - s.spring_force
    nothing
end

"""
Calculate the forces, acting on all particles.
v_wind_tether: out parameter
forces:        out parameter
"""
function inner_loop(s, pos, vel, v_wind_gnd, stiffnes_factor, segments, d_tether)
    for i in 1:length(s.springs)
        p1 = s.springs[i].p1  # First point nr.
        p2 = s.springs[i].p2  # Second point nr.
        height = 0.5 * (pos[p1][3] + pos[p2][3])
        rho = calc_rho(s, height)
        s.v_wind_tether .= calc_wind_factor(s, height) * v_wind_gnd
        calc_particle_forces(s, pos[p1], pos[p2], vel[p1], vel[p2], s.springs[i], stiffnes_factor, segments, d_tether, rho, i)
    end
    nothing
end

"""
pos2, pos3, pos4: position of the kite particles P2, P3, and P4
v2, v3, v4:       velocity of the kite particles P2, P3, and P4
rho:              air density [kg/m^3]
rel_depower:      value between  0.0 and  1.0
rel_steering:     value between -1.0 and +1.0
"""
function calc_aero_forces(s::KPS4, pos, vel, rho, alpha_depower, rel_steering)
    rel_side_area = s.set.rel_side_area/100.0 # defined in percent
    K = 1 - rel_side_area # correction factor for the drag
    pos2, pos3, pos4 = pos[s.set.segments+3], pos[s.set.segments+4], pos[s.set.segments+5]
    v2, v3, v4 = vel[s.set.segments+3], vel[s.set.segments+4], vel[s.set.segments+5]
    va_2, va_3, va_4 = s.v_wind - v2, s.v_wind - v3, s.v_wind - v4
 
    pos_centre = 0.5 * (pos3 + pos4)
    delta = pos2 - pos_centre
    z = -normalize(delta)
    y = normalize(pos3 - pos4)
    x = cross(y, z)
    va_xz2 = va_2 - dot(va_2, y) * y
    va_xy3 = va_3 - dot(va_3, z) * z
    va_xy4 = va_4 - dot(va_4, z) * z

    alpha_2 = (pi - acos(dot(normalize(va_xz2), x)) - alpha_depower) * 360.0 / pi + s.set.alpha_zero
    alpha_3 = (pi - acos(dot(normalize(va_xy3), x)) - rel_steering * KS) * 360.0 / pi + s.set.alpha_ztip
    alpha_4 = (pi - acos(dot(normalize(va_xy4), x)) + rel_steering * KS) * 360.0 / pi + s.set.alpha_ztip

    CL2, CD2 = calc_cl(alpha_2), DRAG_CORR * calc_cd(alpha_2)
    CL3, CD3 = calc_cl(alpha_3), DRAG_CORR * calc_cd(alpha_3)
    CL4, CD4 = calc_cl(alpha_4), DRAG_CORR * calc_cd(alpha_4)

    L2 = (-0.5 * rho * (norm(va_xz2))^2 * s.set.area * CL2) * normalize(cross(va_2, y))
    L3 = (-0.5 * rho * (norm(va_xy3))^2 * s.set.area * rel_side_area * CL3) * normalize(cross(va_3, z))
    L4 = (-0.5 * rho * (norm(va_xy4))^2 * s.set.area * rel_side_area * CL4) * normalize(cross(z, va_4))
    D2 = (-0.5 * K * rho * norm(va_2) * s.set.area * CD2) * va_2
    D3 = (-0.5 * K * rho * norm(va_3) * s.set.area * rel_side_area * CD3) * va_3
    D4 = (-0.5 * K * rho * norm(va_4) * s.set.area * rel_side_area * CD4) * va_4
    s.forces[s.set.segments + 3] .+= (L2 + D2)
    s.forces[s.set.segments + 4] .+= (L3 + D3)
    s.forces[s.set.segments + 5] .+= (L4 + D4)
end

""" 
Calculate the vector res0 using a vector expression, and calculate res1 using a loop
that iterates over all tether segments. 
"""
function loop(s::KPS4, masses, forces, pos, vel, posd, veld, res0, res1)
    L_0      = s.l_tether / s.segments
    # mass_per_meter = s.set.rho_tether * π * (s.set.d_tether/2000.0)^2
    mass_per_meter = 0.011
    res0[1] .= pos[1]
    res1[1] .= vel[1]
    particles = s.set.segments + KITE_PARTICLES + 1
    res0[2:particles] .= vel[2:particles] - posd[2:particles]
    # Compute the masses and forces
    m_tether_particle = mass_per_meter * s.length / L_0
    masses[s.segments+1] .= s.kcu_mass + 0.5 * m_tether_particle
    for i in 1:s.set.segments+1
        masses[i] = m_tether_particle
        s.springs[i].length = s.lenght
        s.springs[i].c_spring = s.c_spring / s.stiffnes_factor
        s.springs[i].damping = s.damping
#     innerLoop2_(pos, vel, vec3[V_wind_gnd], vec3[V_wind_tether], forces, \
#                 scalars[Stiffnes_factor], int(SEGMENTS), D_TETHER)
#     for i in xrange(1, NO_PARTICLES):
#         res1[i] = veld[i] - (G_EARTH - forces[i] / masses[i])
    end
end