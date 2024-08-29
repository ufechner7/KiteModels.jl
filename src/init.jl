const KITE_ANGLE = 3.83 # angle between the kite and the last tether segment due to the mass of the control pod
const PRE_STRESS  = 0.9998   # Multiplier for the initial spring lengths.

# Functions to calculate the inital state vector, the inital masses and initial springs

function init_springs!(s::KPS4)
    l_0     = s.set.l_tether / s.set.segments 
    particles = KiteUtils.get_particles(s.set.height_k, s.set.h_bridle, s.set.width, s.set.m_k)
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

# "get the azimuth of the steering tethers"
# function get_tether_azimuth(width, radius, tip_length, middle_length)
#     α_0 = pi/2 - width/2/radius
#     α_c = α_0 + width*(-2*tip_length + sqrt(2*middle_length^2 + 2*tip_length^2))/(4*(middle_length - tip_length)) / radius
#     distance_c = cos(α_c)*radius
#     α_tether = atan(distance_c/l_tether)
#     return α_tether
# end

# implemented
function get_particles(width, radius, middle_length, tip_length, bridle_center_distance, pos_kite = [ 75., 0., 129.90381057], 
    vec_c=[-15., 0., -25.98076211], v_app=[10.4855, 0, -3.08324])
    # inclination angle of the kite; beta = atan(-pos_kite[2], pos_kite[1]) ???
    beta = pi/2.0

    e_z = normalize(vec_c) # vec_c is the direction of the last two particles
    e_y = normalize(cross(v_app, e_z))
    e_x = normalize(cross(e_y, e_z))

    α_0 = pi/2 - width/2/radius
    α_c = α_0 + width*(-2*tip_length + sqrt(2*middle_length^2 + 2*tip_length^2))/(4*(middle_length - tip_length)) / radius
    α_d = π - α_c

    E = pos_kite
    E_c = pos_kite + e_z * (-bridle_center_distance + radius) # E at center of circle on which the kite shape lies
    C = E_c + e_y*cos(α_c)*radius - e_z*sin(α_c)*radius
    D = E_c + e_y*cos(α_d)*radius - e_z*sin(α_d)*radius

    length(α) = α < π/2 ?
        (tip_length + (middle_length-tip_length)*α*radius/(0.5*width)) :
        (tip_length + (middle_length-tip_length)*(π-α)*radius/(0.5*width))
    P_c = (C+D)./2
    A = P_c - e_x*(length(α_c)*(3/4 - 1/4))

    [E, C, D, A] # important to have the order E = 1, C = 2, D = 3, A = 4
end

# implemented - looks good
function init_springs!(s::KPS4_3L)
    l_0 = s.set.l_tether / s.set.segments
    
    particles = get_particles(s.set.width, s.set.radius, s.set.middle_length, s.set.tip_length, s.set.bridle_center_distance)
    
    # build the tether segments of the three tethers
    k = s.set.e_tether * (s.set.d_tether/2000.0)^2 * pi  / l_0  # Spring stiffness for this spring [N/m]
    c = s.set.damping/l_0                                       # Damping coefficient [Ns/m]
    for i in 1:s.set.segments*3
        s.springs[i] = SP(i, i+3, l_0, k, c)
    end
    
    # build the bridle segments
    for i in 1:KITE_SPRINGS_3L
        p0, p1 = SPRINGS_INPUT_3L[i, 1], SPRINGS_INPUT_3L[i, 2] # bridle points
        l_0 = norm(particles[Int(p1)] - particles[Int(p0)]) * PRE_STRESS
        k = s.set.e_tether * (s.set.d_line/2000.0)^2 * pi / l_0
        p0 += s.num_E-1 # correct the index for the start and end particles of the bridle
        p1 += s.num_E-1
        c = s.set.damping/ l_0
        s.springs[i+s.set.segments*3] = SP(Int(p0), Int(p1), l_0, k, c)
    end
    return s.springs
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

# implemented
function init_masses!(s::KPS4_3L)
    s.masses = zeros(s.num_A)
    l_0 = s.set.l_tether / s.set.segments 
    mass_per_meter = s.set.rho_tether * π * (s.set.d_tether/2000.0)^2
    for i in 4:s.set.segments*3
        s.masses[i]   += l_0 * mass_per_meter
    end
    [s.masses[i] += 0.5 * l_0 * mass_per_meter for i in s.num_E-2:s.num_E]
    s.masses[s.num_E] += 0.5 * s.set.l_bridle * mass_per_meter
    s.masses[s.num_A] += s.set.mass/2
    s.masses[s.num_C] += s.set.mass/4
    s.masses[s.num_D] += s.set.mass/4
    return s.masses 
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
        particles = KiteUtils.get_particles(s.set.height_k, s.set.h_bridle, s.set.width, s.set.m_k)
    else
        particles = KiteUtils.get_particles(s.set.height_k, s.set.h_bridle, s.set.width, s.set.m_k, pos[s.set.segments+1], rotate_in_xz(vec_c, deg2rad(KITE_ANGLE)), s.v_apparent)
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
    for i in eachindex(pos)
        s.pos[i] .= pos[i]
    end
    pos, vel, acc
end

# implemented
function init_pos_vel_acc(s::KPS4_3L, X=zeros(5*s.set.segments+3); delta = 0.0)
    pos = zeros(SVector{s.num_A, KVec3})
    vel = zeros(SVector{s.num_A, KVec3})
    acc = zeros(SVector{s.num_A, KVec3})
    # ground points
    [pos[i] .= [0.0, delta, 0.0] for i in 1:3]

    # middle tether
    sin_el, cos_el = sin(deg2rad(s.set.elevation)), cos(deg2rad(s.set.elevation))
    for (i, j) in enumerate(range(6, step=3, length=s.set.segments))
        radius = i * (s.set.l_tether/s.set.segments)
        pos[j] .= [cos_el*radius + X[i], delta, sin_el*radius + X[s.set.segments+i]]
    end

    # kite points
    vec_c = pos[s.num_E-3] - pos[s.num_E]
    particles = get_particles(s.set.width, s.set.radius, s.set.middle_length, s.set.tip_length, s.set.bridle_center_distance, pos[s.num_E], vec_c, s.v_apparent)
    pos[s.num_A] .= particles[4] + [X[s.set.segments*2+1], 0, X[s.set.segments*2+2]]
    pos[s.num_C] .= particles[2] + X[s.set.segments*2+3 : s.set.segments*2+5]
    pos[s.num_D] .= [pos[s.num_C][1], -pos[s.num_C][2], pos[s.num_C][3]]
    
    # build tether connection points
    e_z = normalize(vec_c)
    distance_c_l = 0.0 # distance between c and left steering line
    # distance_c_l = s.set.tip_length/2 # distance between c and left steering line
    s.tether_lengths[1] = norm(pos[s.num_C] + e_z .* (X[s.set.segments*2+6] + distance_c_l)) # find the right steering tether length
    s.tether_lengths[2] = s.tether_lengths[1]
    pos[s.num_E-2] .= pos[s.num_C] + e_z .* (distance_c_l)
    pos[s.num_E-1] .= pos[s.num_E-2] .* [1.0, -1.0, 1.0]

    # build left and right tether points
    for (i, j) in enumerate(range(4, step=3, length=s.set.segments-1))
        pos[j] .= pos[s.num_E-2] ./ s.set.segments .* i .+
            [X[2*s.set.segments+6+i],
            X[3*s.set.segments+5+i],
            X[4*s.set.segments+4+i]]
        pos[j+1] .= [pos[j][1], -pos[j][2], pos[j][3]]
    end

    # set vel and acc
    for i in 1:s.num_A
        vel[i] .= [delta, delta, delta]
        acc[i] .= [delta, delta, delta]
    end
    
    for i in eachindex(pos)
        s.pos[i] .= pos[i]
    end  

    return pos, vel, acc
end


# Calculate the initial vectors pos and vel. Tether with the initial elevation angle
# se().elevation, particle zero fixed at origin.
# X is a vector of deviations in x and z positions, to be varied to find the inital equilibrium
function init_pos_vel(s::KPS4, X=zeros(2 * (s.set.segments+KITE_PARTICLES)))
    pos, vel, acc = init_pos_vel_acc(s, X; old=true)
    pos, vel
end

function init_pos_vel(s::KPS4_3L, X=zeros(5*s.set.segments+3))
    pos, vel, _ = init_pos_vel_acc(s, X)
    return pos, vel
end

function init_pos(s::KPS4_3L, X=zeros(5*s.set.segments+3); delta=0.0)
    pos_, _, _ = init_pos_vel_acc(s, X; delta=0.0)
    pos = zeros(3, s.num_A)
    [pos[:,i] .= pos_[i] for i in 1:s.num_A]
    return pos
end

function init_inner(s::KPS4, X=zeros(2 * (s.set.segments+KITE_PARTICLES-1)+1); old=false, delta=0.0)
    pos, vel, acc = init_pos_vel_acc(s, X; old=old, delta=delta)
    vcat(pos[2:end], vel[2:end]), vcat(vel[2:end], acc[2:end])
end

function init_inner(s::KPS4_3L, X=zeros(5*s.set.segments+3);delta=0.0)
    pos_, vel_, acc_ = init_pos_vel_acc(s, X; delta=delta)
    # remove last left and right tether point and replace them by the length from C and D
    pos = vcat(
        pos_[4:s.num_E-3], 
        pos_[s.num_E:end],
    )
    len = vcat( # connection length
        norm(pos_[s.num_C]-pos_[s.num_E-2]), # left line connection distance
        norm(pos_[s.num_D]-pos_[s.num_E-1]), # right line connection distance
    )
    vel = vcat(
        vel_[4:s.num_E-3], 
        vel_[s.num_E:end],
    )
    acc = vcat(
        acc_[4:s.num_E-3],
        acc_[s.num_E:end],
    )
    vcat(pos, vel, len, [0,0]), vcat(vel, acc, [0,0], [0,0])
end

# same as above, but returns a tuple of two one dimensional arrays
function init(s::KPS4, X=zeros(2 * (s.set.segments+KITE_PARTICLES-1)+1); old=false, delta=0.0)
    res1_, res2_ = init_inner(s, X; old=old, delta = delta)
    res1, res2  = vcat(reduce(vcat, res1_), [s.l_tether, 0]), vcat(reduce(vcat, res2_),[0,0])
    MVector{6*(s.set.segments+KITE_PARTICLES)+2, SimFloat}(res1), MVector{6*(s.set.segments+KITE_PARTICLES)+2, SimFloat}(res2)
end


# implemented
function init(s::KPS4_3L, X=zeros(5*s.set.segments+3); delta=0.0)
    y_, yd_ = init_inner(s, X; delta = delta)
    y = vcat(reduce(vcat, y_), reduce(vcat,[s.tether_lengths, zeros(3)]))
    yd = vcat(reduce(vcat, yd_), zeros(6))
    MVector{6*(s.num_A-5)+4+6, SimFloat}(y), MVector{6*(s.num_A-5)+4+6, SimFloat}(yd)
end

# rotate a 3d vector around the y axis
function rotate_in_xz(vec, angle)
    result = similar(vec)
    result[1] = cos(angle) * vec[1] - sin(angle) * vec[3]
    result[2] = vec[2]
    result[3] = cos(angle) * vec[3] + sin(angle) * vec[1]
    result
end

