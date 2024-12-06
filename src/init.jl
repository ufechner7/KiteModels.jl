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

function init_springs!(s::KPS4_3L)
    l_0 = s.set.l_tether / s.set.segments
    
    E, C, D, A, _, _ = KiteUtils.get_particles_3l(s.set.width, s.set.radius, 
        s.set.middle_length, s.set.tip_length, s.set.bridle_center_distance)
    particles = [E, C, D, A]
    
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
        p0 += s.num_flap_D # correct the index for the start and end particles of the bridle
        p1 += s.num_flap_D
        c = s.set.damping/ l_0
        s.springs[i+s.set.segments*3] = SP(Int(p0), Int(p1), l_0, k, c)
    end
    return s.springs
end


function init_masses!(s::KPS4)
    MASS_FACTOR = 1.0
    s.masses = zeros(s.set.segments+KITE_PARTICLES+1)
    l_0 = s.set.l_tether / s.set.segments 
    for i in 1:s.set.segments
        s.masses[i]   += 0.5 * l_0 * s.set.rho_tether * (s.set.d_tether/2000.0)^2 * pi
        s.masses[i+1] += 0.5 * l_0 * s.set.rho_tether * (s.set.d_tether/2000.0)^2 * pi
    end
    s.masses[s.set.segments+1] += s.set.kcu_mass * MASS_FACTOR
    k2 = s.set.rel_top_mass * (1.0 - s.set.rel_nose_mass)
    k3 = 0.5 * (1.0 - s.set.rel_top_mass) * (1.0 - s.set.rel_nose_mass)
    k4 = 0.5 * (1.0 - s.set.rel_top_mass) * (1.0 - s.set.rel_nose_mass)
    s.masses[s.set.segments+2] += s.set.rel_nose_mass * s.set.mass * MASS_FACTOR
    s.masses[s.set.segments+3] += k2 * s.set.mass * MASS_FACTOR
    s.masses[s.set.segments+4] += k3 * s.set.mass * MASS_FACTOR
    s.masses[s.set.segments+5] += k4 * s.set.mass * MASS_FACTOR
    s.masses 
end


function init_masses!(s::KPS4_3L)
    s.masses = zeros(s.num_A)
    l_0 = s.set.l_tether / s.set.segments 
    mass_per_meter = s.set.rho_tether * π * (s.set.d_tether/2000.0)^2
    for i in 4:s.set.segments*3
        s.masses[i]   += l_0 * mass_per_meter
    end
    [s.masses[i] += 0.5 * l_0 * mass_per_meter for i in s.num_flap_C:s.num_E]
    s.masses[s.num_E] += 0.5 * s.set.l_bridle * mass_per_meter
    s.masses[s.num_A] += s.set.mass/4
    s.masses[s.num_C] += s.set.mass*3/8
    s.masses[s.num_D] += s.set.mass*3/8
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


function init_pos_vel_acc(s::KPS4_3L; delta = 0.0, α = 45.0)
    pos = zeros(SVector{s.num_A, KVec3})
    vel = zeros(SVector{s.num_A, KVec3})
    acc = zeros(SVector{s.num_A, KVec3})
    # ground points
    [pos[i] .= [0.0, delta, 0.0] for i in 1:3]

    # middle tether
    sin_el, cos_el = sin(deg2rad(s.set.elevation)), cos(deg2rad(s.set.elevation))
    for (i, j) in enumerate(range(6, step=3, length=s.set.segments))
        radius = i * (s.set.l_tether/s.set.segments)
        pos[j] .= [cos_el*radius, delta, sin_el*radius]
    end

    # kite points
    vec_c = pos[s.num_flap_C-1] - pos[s.num_E]
    E, C, D, A, s.α_C, s.kite_length_C = KiteUtils.get_particles_3l(s.set.width, s.set.radius, 
                            s.set.middle_length, s.set.tip_length, s.set.bridle_center_distance, pos[s.num_E], vec_c, s.v_apparent)

    pos[s.num_A] .= A
    pos[s.num_C] .= C
    pos[s.num_D] .= [pos[s.num_C][1], -pos[s.num_C][2], pos[s.num_C][3]]
    
    # build tether connection points
    calc_kite_ref_frame!(s, E, C, D)
    E_C = pos[s.num_E] + s.e_z * (-s.set.bridle_center_distance + s.set.radius) 
    e_r_C = (E_C - pos[s.num_C]) / norm(E_C - pos[s.num_C])
    flap_length = s.kite_length_C/4
    angle_flap_c = 0.0 # distance between c and left steering line
    # distance_c_l = s.set.tip_length/2 # distance between c and left steering line
    pos[s.num_flap_C] .= pos[s.num_C] - s.e_x * flap_length * cos(angle_flap_c) + e_r_C * flap_length * sin(angle_flap_c)
    pos[s.num_flap_D] .= pos[s.num_flap_C] .* [1.0, -1.0, 1.0]
    
    s.tether_lengths[3] = norm(pos[s.num_E])
    s.tether_lengths[1] = 0.0
    # build left and right tether points with degrees of bend α
    l = s.tether_lengths[3]
    @assert s.set.elevation > α > 0.0
    α = deg2rad(α)
    h = l/(2tan(α))
    r = l/(2sin(α))
    for (i, j) in enumerate(range(4, step=3, length=s.set.segments-1))
        # pos[j] .= pos[s.num_flap_C] ./ s.set.segments .* i .+ [(middle_distance)*s.tether_lengths[3]*0.5, 0.0, 0.0]
        γ = -α + 2α*i / s.set.segments
        local_z_minus = l/2 + r * sin(γ)
        local_x = h - r * cos(γ)
        local_y = pos[s.num_flap_C][2] / s.set.segments * i
        pos[j] .= local_z_minus * -s.e_z + local_x * s.e_x + local_y * s.e_y

        s.tether_lengths[1] += norm(pos[j] - pos[j-3])
        pos[j+1] .= [pos[j][1], -pos[j][2], pos[j][3]]
    end
    s.tether_lengths[1] += norm(pos[s.num_flap_C] - pos[s.num_flap_C-3])
    s.tether_lengths[2] = s.tether_lengths[1]

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

function init_pos_vel(s::KPS4_3L; delta=0.0)
    pos, vel, _ = init_pos_vel_acc(s; delta=0.0)
    # pos = Array{Union{Nothing, Float64}}(nothing, 3, s.num_A)
    # vel = Array{Union{Nothing, Float64}}(nothing, 3, s.num_A)
    # [pos[:,i] .= pos_[i] for i in 1:s.num_flap_C-1]
    # [vel[:,i] .= zeros(3) for i in 1:s.num_flap_C-1]
    # [pos[:,i] .= pos_[i] for i in s.num_flap_D+1:s.num_A]
    # [vel[:,i] .= zeros(3) for i in s.num_flap_D+1:s.num_A]
    return pos, vel
end

function init_inner(s::KPS4, X=zeros(2 * (s.set.segments+KITE_PARTICLES-1)+1); old=false, delta=0.0)
    pos, vel, acc = init_pos_vel_acc(s, X; old=old, delta=delta)
    vcat(pos[2:end], vel[2:end]), vcat(vel[2:end], acc[2:end])
end

function init_inner(s::KPS4_3L; delta=0.0)
    pos_, vel_, acc_ = init_pos_vel_acc(s; delta=delta)
    # remove last left and right tether point and replace them by the length from C and D
    pos = vcat(
        pos_[4:s.num_flap_C-1], 
        pos_[s.num_E:end],
    )
    len = vcat( # connection length
        norm(pos_[s.num_C]-pos_[s.num_flap_C]), # left line connection distance
        norm(pos_[s.num_D]-pos_[s.num_flap_D]), # right line connection distance
    )
    vel = vcat(
        vel_[4:s.num_flap_C-1], 
        vel_[s.num_E:end],
    )
    acc = vcat(
        acc_[4:s.num_flap_C-1],
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

function init(s::KPS4_3L; delta=0.0)
    y_, yd_ = init_inner(s; delta = delta)
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

# rotate a 3d vector around the z axis
function rotate_in_yx(vec, angle)
    result = similar(vec)
    result[1] = cos(angle) * vec[1] + sin(angle) * vec[2]
    result[2] = cos(angle) * vec[2] - sin(angle) * vec[1]
    result[3] = vec[3]
    result
end