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
                       1.    2.   -1. # s1, p7, p8
                       2.    3.   -1. # s2, p8, p9
                       3.    4.   -1. # s3, p9, p10
                       3.    5.   -1. # s4, p9, p11
                       4.    1.   -1. # s5, p10, p7
                       5.    1.   -1. # s6, p11, p7 not in diagram
                       4.    5.   -1. # s7, p10, p11
                       4.    2.   -1. # s9, p10, p8 
                       5.    2.   -1.]# s9, p11, p8

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
const KITE_ANGLE = 3.8 # angle between the kite and the last tether segment due to the mass of the control pod
const DELTA_MAX = 30.0
const USE_NOMAD = true
const MAX_INTER  = 100  # max interations for steady state finder
const PRE_STRESS  = 0.9998   # Multiplier for the initial spring lengths.
const KS = deg2rad(16.565 * 1.064 * 0.875 * 1.033 * 0.9757 * 1.083)  # max steering
const DRAG_CORR = 0.93       # correction of the drag for the 4-point model
# const X0 = [0.28145470281885937, 0.23227171443921474, -0.14839155793746792, -0.8611709892573215, -1.9064957057144336, -3.2593647930763865, -0.07481136068131161, -0.0320181798695833,  0.12347499137210405,  0.386732547359143,   0.7527967363707873, 1.2094642620268343]
# const X00 = [1.2920877908591142, 1.805784436840238, 1.9811991362643417, 2.1624168969962803, 2.308053915351016, 2.426323438597114, -0.347399490607679, -0.3382081199552657, -0.47968612691009005, -0.48986133189890857, -0.6729736083472052, -0.7322501981033889, -0.7935456484742676, -0.8423258191651944, -0.8811464437359281, 0.04970303818984063, 0.10393814212523027, 0.16883269717897145, -0.003115768195860734]
# const X00 = [1.286733,  1.807743,   2.019513,   2.172045,   2.305881,   2.665952,  -0.394539,  -0.267354,  -0.746302,  -0.487799,  -0.673806,  -0.7461,    -0.797136,  -0.84151,   -0.967774,  -0.012006,   0.066886,   0.228953,  -0.038472]
#const X00 = [1.29044,    1.805875,   2.024852,   2.187221,   2.31497,    2.669564,  -0.393218,  -0.265462,  -0.757762,  -0.487393,  -0.669517,  -0.742705,  -0.795564,  -0.836682,  -0.960051,  -0.013596,   0.070221,   0.233618,  -0.038532]
#const X00 = [1.310236,   1.826226,   2.058625,   2.224326,   2.324677,   2.801553,  -0.400658,  -0.275782,  -0.767225,  -0.49256,   -0.671764,  -0.746666, -0.797588,  -0.825812,  -0.991318,  -0.016655,   0.070457,   0.226817,  -0.037857 ]
#const X00 = [1.9616380996056078, 3.0831297183579403, 3.923474337317383, 4.709346372942605, 5.456808834543749, 6.635620447165995, -0.5304744112903841, -0.3953419401018206, -0.9395790446300472, -0.7700591568439518, -1.1820423830524633, -1.4821917813075025, -1.7611277971401813, -2.0251961180248923, -2.4604625120224237, 0.0014807397938664249, 0.11049910426445765, 0.2876133948013786, -0.0384747365700002]
#const X00 = [1.0989198969424006, 1.6148835475875825, 1.8901129999579414, 2.0875735022084823, 2.2261028201663104, 2.582785078998004, -0.24569831479342286, -0.11401194093941638, -0.6258988098871039, -0.3676736927862901, -0.5088013058266545, -0.5612169010415446, -0.5855688217128946, -0.5888548982433978, -0.6709683952476885, -0.0354902538625784, 0.048556577827989514, 0.2004032421566736, -0.04191614720724987]
const X00 = [0.894146,   1.528959,   1.832319,   2.133604,   2.213243,   2.096041,  -0.242634,  -0.121017,  -0.586176,  -0.300887,  -0.501362,  -0.578397,  -0.654672,  -0.651483,  -0.579665,  -0.034478,   0.037224,   0.17459,   -0.039631]
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
    "drag force of kite and bridle; output of calc_aero_forces"
    drag_force::T =       zeros(S, 3)
    "lift force of the kite; output of calc_aero_forces"
    lift_force::T =       zeros(S, 3)    
    "spring force of the current tether segment, output of calc_particle_forces"
    spring_force::T =     zeros(S, 3)
    segment::T =          zeros(S, 3)
    v_kite::T =           zeros(S, 3)        
    res1::SVector{P, KVec3} = zeros(SVector{P, KVec3})
    res2::SVector{P, KVec3} = zeros(SVector{P, KVec3})
    pos::SVector{P, KVec3} = zeros(SVector{P, KVec3})
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
    v_reel_out::S =        0.0
    last_v_reel_out::S =   0.0
    l_tether::S =          0.0
    rho::S =               0.0
    depower::S =           0.0
    steering::S =          0.0
    stiffness_factor::S =  1.0
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
    s.segment_length = s.l_tether / s.set.segments
    s.v_kite = zeros(SimFloat, 3)
    s.beta = deg2rad(s.set.elevation)
    init_masses(s)
    init_springs(s)
    for i in 1:se().segments + KiteModels.KITE_PARTICLES + 1 
        s.forces[i] .= zeros(3)
    end
    s.drag_force .= [0.0, 0, 0]
    s.lift_force .= [0.0, 0, 0]
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
Parameters:
- mk: relative nose distance
"""
function get_particles(height_k, height_b, width, m_k, pos_pod= [ 75., 0., 129.90381057], vec_c=[-15., 0., -25.98076211], v_app=[10.4855, 0, -3.08324])
    # inclination angle of the kite; beta = atan(-pos_kite[2], pos_kite[1]) ???
    beta = pi/2.0
    x, y, z = initial_kite_ref_frame(vec_c, v_app)

    h_kx = height_k * cos(beta); # print 'h_kx: ', h_kx
    h_kz = height_k * sin(beta); # print 'h_kz: ', h_kz
    h_bx = height_b * cos(beta)
    h_bz = height_b * sin(beta)
    pos_kite = pos_pod - (h_kz + h_bz) * z + (h_kx + h_bx) * x  # top,        poing B in diagram
    pos3 = pos_kite + h_kz * z + 0.5 * width * y + h_kx * x     # side point, point C in diagram
    pos1 = pos_kite + h_kz * z + (h_kx + width * m_k) * x       # nose,       point A in diagram
    pos4 = pos_kite + h_kz * z - 0.5 * width * y + h_kx * x     # side point, point D in diagram
    pos0 = pos_kite + (h_kz + h_bz) * z + (h_kx + h_bx) * x     # equal to pos_pod, P_KCU in diagram
    # println("norm(pos1-pos3) $(norm(pos1-pos3))")             # S2 p8 p9 
    # println("norm(pos4-pos3) $(norm(pos4-pos3))")             # S3 p9 p10
    # println("norm(pos1-pos4) $(norm(pos1-pos4))")             # S9 p11 p8
    [zeros(3), pos0, pos1, pos3, pos4, pos_kite] # 0, p7, p8, p9, p10, p11
end

function calc_height(s::KPS4)
    pos_kite = 0.5 * (s.pos[s.set.segments+4] + s.pos[s.set.segments+5])
    pos_kite[3]
end

function init_springs(s::KPS4)
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
                # if j+s.set.segments-1 == 15
                #     println("s9: $(l_0), $p0, $p1")
                # end
                s.springs[j+s.set.segments-1] = SP(Int(p0), Int(p1), l_0, k, c)
            end
        end
    end
    s.springs
end

function init_masses(s::KPS4)
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

function init_pos_vel_acc(s::KPS4, X=zeros(2 * (s.set.segments+KITE_PARTICLES-1)+1); old=false)
    delta = 1e-6
    pos = zeros(SVector{s.set.segments+1+KITE_PARTICLES, KVec3})
    vel = zeros(SVector{s.set.segments+1+KITE_PARTICLES, KVec3})
    acc = zeros(SVector{s.set.segments+1+KITE_PARTICLES, KVec3})
    pos[1] .= [0.0, delta, 0.0]
    vel[1] .= [delta, delta, delta]
    acc[1] .= [delta, delta, delta]
    sin_el, cos_el = sin(s.set.elevation / 180.0 * pi), cos(s.set.elevation / 180.0 * pi)
    for i in 1:s.set.segments
        radius = -i * (s.set.l_tether/s.set.segments)
        if old
            pos[i+1] .= [-cos_el * radius + X[i], delta, -sin_el * radius + X[s.set.segments+KITE_PARTICLES-1+i]]
        else
            # pos[i+1] .= [-cos_el * radius + X[i] + X0[i], delta, -sin_el * radius + X[s.set.segments+KITE_PARTICLES-1+i] + X0[s.set.segments+i]]
            pos[i+1] .= [-cos_el * radius + X[i], delta, -sin_el * radius + X[s.set.segments+KITE_PARTICLES-1+i]]
        end
        vel[i+1] .= [delta, delta, 0]
        acc[i+1] .= [delta, delta, -9.81]
    end
    vec_c = pos[6] - pos[7]
    if old
        particles = get_particles(s.set.height_k, s.set.h_bridle, s.set.width, s.set.m_k)
    else
        particles = get_particles(s.set.height_k, s.set.h_bridle, s.set.width, s.set.m_k, pos[s.set.segments+1], rotate_in_xz(vec_c, deg2rad(KITE_ANGLE)), s.v_apparent)
    end
    j = 1
    for i in [1,2,4] # set p8, p9, p11
        pos[s.set.segments+1+i] .= particles[i+2] + [X[s.set.segments+j], 0, X[2*s.set.segments+KITE_PARTICLES-1+j]]
        vel[s.set.segments+1+i] .= [delta, delta, delta]
        acc[s.set.segments+1+i] .= [delta, delta, -9.81]
        j +=1
    end
    # set p9=C and p10=D
    # x and z component of the right and left particle must be equal
    pos[s.set.segments+1+3][1] = pos[s.set.segments+1+2][1]  # D.x = C.x
    pos[s.set.segments+1+3][3] = pos[s.set.segments+1+2][3]  # D.z = C.z
    pos[s.set.segments+1+2][2] += X[end]                     # Y position of point C
    pos[s.set.segments+1+3][2] = -pos[s.set.segments+1+2][2] # Y position of point D
    for i in 1:length(pos)
        s.pos[i] .= pos[i]
    end
    # println("pos[8], pos[10]: $(pos[8]), $(pos[10])")
    pos, vel, acc
end

function init(s::KPS4, X=zeros(2 * (s.set.segments+KITE_PARTICLES-1)+1); old=false)
    pos, vel, acc = init_pos_vel_acc(s, X; old=old)
    vcat(pos, vel), vcat(vel, acc)
end

# same as above, but returns a tuple of two one dimensional arrays
function init_flat(s::KPS4, X=zeros(2 * (s.set.segments+KITE_PARTICLES-1)); old=false)
    res1_, res2_ = init(s, X; old=old)
    res1, res2  = reduce(vcat, res1_), reduce(vcat, res2_)
    # append the initial reel-out length and it's derivative
    res1 = vcat(res1, SVector(s.set.l_tether, s.set.v_reel_out))
    res2 = vcat(res2, SVector(s.set.v_reel_out, 1e-6))
    res1, res2
end

""" 
Calculate the drag force of the tether segment, defined by the parameters pos1, pos2, vel1 and vel2
and distribute it equally on the two particles, that are attached to the segment.
The result is stored in the array s.forces. 
"""
@inline function calc_particle_forces(s, pos1, pos2, vel1, vel2, spring, segments, d_tether, rho, i)
    l_0 = spring.length # Unstressed length
    k = spring.c_spring * s.stiffness_factor  # Spring constant
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
@inline function inner_loop(s, pos, vel, v_wind_gnd, segments, d_tether)
    for i in 1:length(s.springs)
        p1 = s.springs[i].p1  # First point nr.
        p2 = s.springs[i].p2  # Second point nr.
        height = 0.5 * (pos[p1][3] + pos[p2][3])
        rho = calc_rho(s, height)
        s.v_wind_tether .= calc_wind_factor(s, height) * v_wind_gnd
        calc_particle_forces(s, pos[p1], pos[p2], vel[p1], vel[p2], s.springs[i], segments, d_tether, rho, i)
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

    alpha_2 = (pi - acos(max(dot(normalize(va_xz2), x), -1.0)) - alpha_depower) * 180.0 / pi + s.set.alpha_zero
    alpha_3 = (pi - acos(max(dot(normalize(va_xy3), x), -1.0)) - rel_steering * KS) * 180.0 / pi + s.set.alpha_ztip
    alpha_4 = (pi - acos(max(dot(normalize(va_xy4), x), -1.0)) + rel_steering * KS) * 180.0 / pi + s.set.alpha_ztip

    CL2, CD2 = calc_cl(alpha_2), DRAG_CORR * calc_cd(alpha_2)
    CL3, CD3 = calc_cl(alpha_3), DRAG_CORR * calc_cd(alpha_3)
    CL4, CD4 = calc_cl(alpha_4), DRAG_CORR * calc_cd(alpha_4)
    L2 = (-0.5 * rho * (norm(va_xz2))^2 * s.set.area * CL2) * normalize(cross(va_2, y))
    L3 = (-0.5 * rho * (norm(va_xy3))^2 * s.set.area * rel_side_area * CL3) * normalize(cross(va_3, z))
    L4 = (-0.5 * rho * (norm(va_xy4))^2 * s.set.area * rel_side_area * CL4) * normalize(cross(z, va_4))
    D2 = (-0.5 * K * rho * norm(va_2) * s.set.area * CD2) * va_2
    D3 = (-0.5 * K * rho * norm(va_3) * s.set.area * rel_side_area * CD3) * va_3
    D4 = (-0.5 * K * rho * norm(va_4) * s.set.area * rel_side_area * CD4) * va_4
    s.lift_force .= L2
    # println("L3, L4: $(norm(L3)), $(norm(L4))")
    # println("D2, D3, D4: $D2, $D3, $D4")
    s.drag_force .= D2 + D3 + D4
    s.forces[s.set.segments + 3] .+= (L2 + D2)
    s.forces[s.set.segments + 4] .+= (L3 + D3)
    s.forces[s.set.segments + 5] .+= (L4 + D4)
end

""" 
Calculate the vector res1 and calculate res2 using loops
that iterate over all tether segments. 
"""
function loop(s::KPS4, pos, vel, posd, veld)
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
    inner_loop(s, pos, vel, s.v_wind_gnd, s.set.segments, s.set.d_tether/1000.0)
    for i in 2:particles
        s.res2[i] .= veld[i] - (SVector(0, 0, -G_EARTH) - s.forces[i] / s.masses[i])
    end
end

"""
    residual!(res, yd, y::MVector{S, SimFloat}, s::KPS3, time) where S

    N-point tether model, one point kite at the top:
    Inputs:
    State vector y   = pos1,  pos2, ... , posn,  vel1,  vel2, . .., veln,  length, v_reel_out
    Derivative   yd  = posd1, posd2, ..., posdn, veld1, veld2, ..., veldn, lengthd, v_reel_outd
    Output:
    Residual     res = res1, res2 = vel1-posd1,  ..., veld1-acc1, ...

    Additional parameters:
    s: Struct with work variables, type KPS3
    S: The dimension of the state vector
The number of the point masses of the model N = (S-2)/6, the state of each point 
is represented by two 3 element vectors.
"""
function residual!(res, yd, y::MVector{S, SimFloat}, s::KPS4, time) where S
    T = S-2 # T: three times the number of particles including the origin
    segments = div(T,6) - KITE_PARTICLES - 1
    
    # Reset the force vector to zero.
    for i in 1:segments + KITE_PARTICLES + 1 
        s.forces[i] .= SVector(0.0, 0, 0)
    end
    length = y[end-1]
    v_reel_out = y[end]
    lengthd = yd[end-1]
    v_reel_outd = yd[end]

    # unpack the vectors y and yd
    ys  = @view y[1:T]
    yds = @view yd[1:T]
    part  = reshape(SVector{T}(ys),  Size(3, div(T,6), 2))
    partd = reshape(SVector{T}(yds), Size(3, div(T,6), 2))
    pos1, vel1 = part[:,:,1], part[:,:,2]
    pos = SVector{div(T,6)}(SVector(pos1[:,i]) for i in 1:div(T,6))
    vel = SVector{div(T,6)}(SVector(vel1[:,i]) for i in 1:div(T,6))
    posd1, veld1 = partd[:,:,1], partd[:,:,2]
    posd = SVector{div(T,6)}(SVector(posd1[:,i]) for i in 1:div(T,6))
    veld = SVector{div(T,6)}(SVector(veld1[:,i]) for i in 1:div(T,6))

    # core calculations
    calc_aero_forces(s, pos, vel, s.rho, s.alpha_depower, s.steering)
    loop(s, pos, vel, posd, veld)

    # copy and flatten result
    for i in 2:div(T,6)
        for j in 1:3
            res[3*(i-2)+j]                = s.res1[i][j]
            res[3*(div(T,6)-1)+3*(i-2)+j] = s.res2[i][j]
        end
    end
    if norm(res) < 1e5
        for i in 1:div(T,6)
            @inbounds s.pos[i] .= pos[i]
        end
    end

    # winch not yet integrated
    res[end-1] = 0.0
    res[end]   = 0.0
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
            s.spring_force .= k *  (norm1 - l_0) 
        else 
            s.spring_force .= k1 *  (norm1 - l_0)
        end
        forces[i+s.set.segments] = norm(s.spring_force)
        if norm(s.spring_force) > 4000.0
            println("Bridle brakes for spring $i !")
        end
    end
    forces
end

"""
    find_steady_state(s::KPS4, prn=false)

Find an initial equilibrium, based on the inital parameters
`l_tether`, elevation and `v_reel_out`.
"""
function find_steady_state(s::KPS4, prn=false)
    res = zeros(MVector{6*(s.set.segments+KITE_PARTICLES)+2, SimFloat})
    n_ones = ones(SimFloat, 2*(s.set.segments+KITE_PARTICLES-1)+1)
    iter = 0

    # helper function for the steady state finder
    function test_initial_condition!(F, x::Vector)
        x1 = copy(x)
        y0, yd0 = init_flat(s, x1)
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
        # copy the acceleration of point C in y direction
        i = s.set.segments+3
        F[end]                                 = res[2 + 3*(i-1) + 3*(s.set.segments+KITE_PARTICLES)] 
        iter += 1
        return nothing 
    end
    function f(x)
      F = zeros(SimFloat, 2*(s.set.segments+KITE_PARTICLES-1)+1)
      test_initial_condition!(F, x)
      # reduce the weight of the tether particles
      for i in 1:s.set.segments       
          F[i] *= 0.25 / (i/(2*s.set.segments))
          F[i+KITE_PARTICLES-1] *= 0.5 / (i/(2*s.set.segments))
      end    
      F[end] *= 5.0  
      F
    end
    function eval_fct(x)
        bb_outputs = [norm(f(x))]
        success = ! isnan(bb_outputs[1]) # && all(KiteModels.spring_forces(s) .> 0)
        count_eval = true
        success, count_eval, bb_outputs
    end
    pb = NomadProblem(2*(s.set.segments+KITE_PARTICLES-1)+1, # number of inputs of the blackbox
                    1, # number of outputs of the blackbox
                    ["OBJ"], # type of outputs of the blackbox
                    eval_fct;
                    lower_bound = -DELTA_MAX * n_ones,
                    upper_bound =  DELTA_MAX * n_ones)
    pb.options.max_bb_eval = MAX_INTER * 10
    if prn println("\nStarted function test_nlsolve...") end
    if USE_NOMAD
        # result = solve(pb, zeros(SimFloat, 2*(s.set.segments+KITE_PARTICLES-1)+1))
        result = solve(pb, X00)       
        init(s, result[1])
    else
        results = nlsolve(test_initial_condition!, X00, autoscale=true, iterations=MAX_INTER)
        if prn println("\nresult: $results") end
        init(s, results.zero)
    end
end

# rotate a 3d vector around the y axis
function rotate_in_xz(vec, angle)
    result = similar(vec)
    result[1] = cos(angle) * vec[1] - sin(angle) * vec[3]
    result[2] = vec[2]
    result[3] = cos(angle) * vec[3] + sin(angle) * vec[1]
    result
end
