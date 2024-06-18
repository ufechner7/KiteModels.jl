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

# implemented
function get_particles(height_k, height_b, width, m_k, pos_E = [ 75., 0., 129.90381057], 
    vec_c=[-15., 0., -25.98076211], v_app=[10.4855, 0, -3.08324])
    # inclination angle of the kite; beta = atan(-pos_kite[2], pos_kite[1]) ???
    beta = pi/2.0

    e_z = normalize(vec_c) # vec_c is the direction of the last two particles
    e_y = normalize(cross(v_app, vec_c))
    e_x = normalize(cross(y, vec_c))

    α_c = α_0 + s.set.width*(-2*s.set.tip_length + sqrt(2*s.set.middle_length^2 + 2*s.set*tip_length^2))/(4*(s.set.middle_length - s.set.tip_length)) / s.set.radius
    α_d = π - α_c

    E = pos_E
    C = E + e_y*cos(α_c)*r - e_z*sin(α_c)*r
    D = E + e_y*cos(α_d)*r - e_z*sin(α_d)*r

    length(α) = α < π/2 ?
        (s.set.tip_length + (s.set.middle_length-s.set.tip_length)*α*s.set.radius/(0.5*s.set.width)) :
        (s.set.tip_length + (s.set.middle_length-s.set.tip_length)*(π-α)*s.set.radius/(0.5*s.set.width))
    P_c = (C+D)./2
    A = P_c + e_x.*length(α_c)*(-3/4 + 1/4)

    [zeros(3), E, C, D, A] # 0, E, A, C, D
end

# implemented
function init_springs!(s::KPS4_3L)
    l_0 = s.set.l_tether / s.set.segments
    
    particles = get_particles(s.set.height_k, s.set.h_bridle, s.set.width, s.set.m_k)
    
    k = s.set.e_tether * (s.set.d_tether/2000.0)^2 * pi  / l_0  # Spring stiffness for this spring [N/m]
    c = s.set.damping/l_0                                       # Damping coefficient [Ns/m]

    # build the tether segments of the three tethers
    [s.springs[i] = SP(1, i, l_0, k, c) for i in 2:4]
    for i in 2:3:s.set.segments*3
        s.springs[i+1] = SP(i, i+3, l_0, k, c) # left line
        s.springs[i+2] = SP(i+1, i+4, l_0, k, c) # right line
        s.springs[i] = SP(i+2, i+5, l_0, k, c) # middle line
    end

    # build the bridle segments
    for j in 4:size(SPRINGS_INPUT_3L)[1] # 4th input is the first bridle segment
        p0, p1 = SPRINGS_INPUT_3L[j, 1]+1, SPRINGS_INPUT_3L[j, 2]+1 # point 0 and 1
        if SPRINGS_INPUT_3L[j, 3] == -1
            l_0 = norm(particles[Int(p1)] - particles[Int(p0)]) * PRE_STRESS
            k = s.set.e_tether * (s.set.d_line/2000.0)^2 * pi / l_0
            p0 += s.set.segments*3 - 1 # correct the index for the start and end particles of the bridle
            p1 += s.set.segments*3 - 1
            c = s.set.damping/ l_0
            s.springs[j+s.set.segments*3-1] = SP(Int(p0), Int(p1), l_0, k, c)
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

# implemented
function init_masses!(s::KPS4_3L)
    s.masses = zeros(s.set.segments+KITE_PARTICLES_3L+1)
    l_0 = s.set.l_tether / s.set.segments 
    for i in 1:s.set.segments
        s.masses[i]   += 0.5 * l_0 * s.set.rho_tether * (s.set.d_tether/2000.0)^2 * pi
        s.masses[i+1] += 0.5 * l_0 * s.set.rho_tether * (s.set.d_tether/2000.0)^2 * pi
    end
    m_a = s.set.mass/2
    m_c = s.set.mass/4
    m_d = s.set.mass/4
    s.masses[s.set.segments+4] += m_a
    s.masses[s.set.segments+2] += m_c
    s.masses[s.set.segments+3] += m_d
    s.masses 
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
function init_pos_vel_acc(s::KPS4_3L, X=zeros(2 * (s.set.segments+KITE_PARTICLES_3L)+1); delta = 0.0)
    pos = zeros(SVector{s.set.segments*3+KITE_PARTICLES_3L, KVec3})
    vel = zeros(SVector{s.set.segments*3+KITE_PARTICLES_3L, KVec3})
    acc = zeros(SVector{s.set.segments*3+KITE_PARTICLES_3L, KVec3})
    pos[1] .= [0.0, delta, 0.0]
    vel[1] .= [delta, delta, delta]
    acc[1] .= [delta, delta, delta]
    sin_el, cos_el = sin(deg2rad(s.set.elevation)), cos(deg2rad(s.set.elevation))

    radius = -1 * (s.set.l_tether/s.set.segments)
    pos[2:4] .= [[-cos_el * radius + X[i], delta, -sin_el * radius + X[s.set.segments*3+KITE_PARTICLES_3L-1+i]]]
    vel[2:4] .= [[delta, delta, 0]]
    acc[2:4] .= [[delta, delta, -9.81]]
    for i in 2:s.set.segments*3-2
        radius = -i * (s.set.l_tether/s.set.segments)
        pos[i+3] .= [-cos_el * radius + X[i], delta, -sin_el * radius + X[s.set.segments*3+KITE_PARTICLES_3L-1+i]]
        vel[i+3] .= [delta, delta, 0]
        acc[i+3] .= [delta, delta, -9.81]
    end

    vec_c = pos[s.set.segments*3-2] - pos[s.set.segments*3+1]
    particles = get_particles(s.set.height_k, s.set.h_bridle, s.set.width, s.set.m_k, pos[s.set.segments*3+1], rotate_in_xz(vec_c, deg2rad(KITE_ANGLE_3L)), s.v_apparent)

    for i in [2,4] # set A, B, C
        pos[s.set.segments*3+i] .= particles[i+1] + [X[s.set.segments*3+i-1], 0, X[2*s.set.segments*3+KITE_PARTICLES_3L-2+i]]
        vel[s.set.segments*3+i] .= [delta, delta, delta]
        acc[s.set.segments*3+i] .= [delta, delta, -9.81]
    end
    acc[s.set.segments*3+3] .= [delta, delta, -9.81]
    vel[s.set.segments*3+3] .= [delta, delta, delta]
    # set p10=C and p11=D
    # x and z component of the right and left particle must be equal
    pos[s.set.segments*3+2][1] = pos[s.set.segments*3+1+3][1]  # D.x = C.x
    pos[s.set.segments*3+3][3] = pos[s.set.segments*3+1+3][3]  # D.z = C.z
    pos[s.set.segments*3+2][2] += X[end]                     # Y position of point C
    pos[s.set.segments*3+3][2] = -pos[s.set.segments*3+1+3][2] # Y position of point D
    for i in eachindex(pos)
        s.pos[i] .= pos[i]
    end
    pos, vel, acc
end


# Calculate the initial vectors pos and vel. Tether with the initial elevation angle
# se().elevation, particle zero fixed at origin.
# X is a vector of deviations in x and z positions, to be varied to find the inital equilibrium
function init_pos_vel(s::KPS4, X=zeros(2 * (s.set.segments+KITE_PARTICLES)))
    pos, vel, acc = init_pos_vel_acc(s, X; old=true)
    pos, vel
end

function init_inner(s::KPS4, X=zeros(2 * (s.set.segments+KITE_PARTICLES-1)+1); old=false, delta=0.0)
    pos, vel, acc = init_pos_vel_acc(s, X; old=old, delta=delta)
    vcat(pos[2:end], vel[2:end]), vcat(vel[2:end], acc[2:end])
end

function init_inner(s::KPS4_3L, X=zeros(2 * (s.set.segments+KITE_PARTICLES_3L-1)+1);delta=0.0)
    pos, vel, acc = init_pos_vel_acc(s, X; delta=delta)
    vcat(pos[2:end], vel[2:end]), vcat(vel[2:end], acc[2:end])
end

# same as above, but returns a tuple of two one dimensional arrays
function init(s::KPS4, X=zeros(2 * (s.set.segments+KITE_PARTICLES-1)+1); old=false, delta=0.0)
    res1_, res2_ = init_inner(s, X; old=old, delta = delta)
    res1, res2  = vcat(reduce(vcat, res1_), [s.l_tether, 0]), vcat(reduce(vcat, res2_),[0,0])
    MVector{6*(s.set.segments+KITE_PARTICLES)+2, SimFloat}(res1), MVector{6*(s.set.segments+KITE_PARTICLES)+2, SimFloat}(res2)
end


"""
State vector y   = pos1,  pos2, ... , posn,  vel1,  vel2, . .., veln,  length1, length2, length3, v_reel_out1, v_reel_out2, v_reel_out3
Derivative   yd  = posd1, posd2, ..., posdn, veld1, veld2, ..., veldn, lengthd1, lengthd2, lengthd3, v_reel_outd1, v_reel_outd2, v_reel_outd3
"""
# same as above, but returns a tuple of two one dimensional arrays
function init(s::KPS4_3L, X=zeros(2 * (s.set.segments+KITE_PARTICLES_3L-1)+1); delta=0.0)
    y_, yd_ = init_inner(s, X; delta = delta)
    y, yd  = vcat(reduce(vcat, y_), [s.l_tether, 0]), vcat(reduce(vcat, yd_),[0,0])
    MVector{6*(s.set.segments+KITE_PARTICLES_3L)+2, SimFloat}(y), MVector{6*(s.set.segments+KITE_PARTICLES_3L)+2, SimFloat}(yd)
end

# rotate a 3d vector around the y axis
function rotate_in_xz(vec, angle)
    result = similar(vec)
    result[1] = cos(angle) * vec[1] - sin(angle) * vec[3]
    result[2] = vec[2]
    result[3] = cos(angle) * vec[3] + sin(angle) * vec[1]
    result
end