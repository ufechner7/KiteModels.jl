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

function init_springs!(s::KPS4_3L)
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

function init_masses!(s::KPS4_3L)
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

function init_pos_vel_acc(s::KPS4_3L, X=zeros(2 * (s.set.segments+KITE_PARTICLES)+1); old=false, delta = 0.0)
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

function init_inner(s::KPS4_3L, X=zeros(2 * (s.set.segments+KITE_PARTICLES-1)+1); old=false, delta=0.0)
    pos, vel, acc = init_pos_vel_acc(s, X; old=old, delta=delta)
    vcat(pos[2:end], vel[2:end]), vcat(vel[2:end], acc[2:end])
end

# same as above, but returns a tuple of two one dimensional arrays
function init(s::KPS4, X=zeros(2 * (s.set.segments+KITE_PARTICLES-1)+1); old=false, delta=0.0)
    res1_, res2_ = init_inner(s, X; old=old, delta = delta)
    res1, res2  = vcat(reduce(vcat, res1_), [s.l_tether, 0]), vcat(reduce(vcat, res2_),[0,0])
    MVector{6*(s.set.segments+KITE_PARTICLES)+2, SimFloat}(res1), MVector{6*(s.set.segments+KITE_PARTICLES)+2, SimFloat}(res2)
end

# same as above, but returns a tuple of two one dimensional arrays
function init(s::KPS4_3L, X=zeros(2 * (s.set.segments+KITE_PARTICLES-1)+1); old=false, delta=0.0)
    res1_, res2_ = init_inner(s, X; old=old, delta = delta)
    res1, res2  = vcat(reduce(vcat, res1_), [s.l_tether, 0]), vcat(reduce(vcat, res2_),[0,0])
    MVector{6*(s.set.segments+KITE_PARTICLES)+2, SimFloat}(res1), MVector{6*(s.set.segments+KITE_PARTICLES)+2, SimFloat}(res2)
end

# rotate a 3d vector around the y axis
function rotate_in_xz(vec, angle)
    result = similar(vec)
    result[1] = cos(angle) * vec[1] - sin(angle) * vec[3]
    result[2] = vec[2]
    result[3] = cos(angle) * vec[3] + sin(angle) * vec[1]
    result
end