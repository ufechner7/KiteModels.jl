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

# function init_springs!(s::KPSQ)
#     l_0 = s.set.l_tether / s.set.segments
    
#     E, C, D, A, _, _ = KiteUtils.get_particles_3l(s.set.width, s.set.radius, 
#         s.set.middle_length, s.set.tip_length, s.set.bridle_center_distance)
#     particles = [E, C, D, A]
    
#     # build the tether segments of the three tethers
#     k = s.set.e_tether * (s.set.d_tether/2000.0)^2 * pi  / l_0  # Spring stiffness for this spring [N/m]
#     c = s.set.damping/l_0                                       # Damping coefficient [Ns/m]
#     for i in 1:s.set.segments*3
#         s.springs[i] = SP(i, i+3, l_0, k, c)
#     end
#     return s.springs
# end


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


function init_masses!(s::KPSQ)
    s.masses = zeros(s.i_C)
    l_0 = s.set.l_tether / s.set.segments 
    mass_per_meter = s.set.rho_tether * π * (s.set.d_tether/2000.0)^2
    for i in 4:s.i_C
        s.masses[i]   += l_0 * mass_per_meter
    end
    [s.masses[i] += 0.5 * l_0 * mass_per_meter for i in s.i_A:s.i_C]
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
        particles = KiteUtils.get_particles(s.set.height_k, s.set.h_bridle, s.set.width, s.set.m_k, 
                              pos[s.set.segments+1], rotate_around_y(vec_c, deg2rad(KITE_ANGLE)), s.v_apparent)
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

function calc_inertia!(s::KPSQ)
    segs = 100
    mass_per_area = s.set.mass / ((s.set.middle_length + s.set.tip_length) * 0.5 * s.set.width)
    
    # First pass - calculate COM
    total_mass = 0.0
    s.pos_circle_center_b .= 0.0 # translation from kite COM to kite circle center
    pos = zeros(3, 2segs)
    mass = zeros(2segs)
    
    @assert s.α_l != 0.0
    α_middle    = π/2
    dα          = (α_middle - s.α_l) / segs
    for i in 1:2segs
        if i <= segs
            α = s.α_l + -dα/2 + i * dα
            kite_length = s.set.tip_length + (s.set.middle_length-s.set.tip_length) * (α - s.α_l) / (π/2 - s.α_l)
        else
            α = pi - (s.α_l + -dα/2 + (i-segs) * dα)
            kite_length = s.set.tip_length + (s.set.middle_length-s.set.tip_length) * (π - s.α_l - α) / (π/2 - s.α_l)
        end
        
        area = kite_length * s.set.width/(2segs)
        mass[i] = area * mass_per_area
        pos[:, i] = [-0.5 * kite_length, cos(α) * s.set.radius, sin(α) * s.set.radius]
        s.pos_circle_center_b .-= mass[i] * pos[:, i]
    end
    s.pos_circle_center_b ./= sum(mass)
    
    # Calculate full inertia tensor relative to COM
    I = zeros(3,3)
    for i in eachindex(mass)
        pos[:, i] -= s.pos_circle_center_b
        r = @views pos[:, i]
        m = mass[i]
        
        # Diagonal terms
        I[1,1] += m * (r[2]^2 + r[3]^2)  # Ixx
        I[2,2] += m * (r[1]^2 + r[3]^2)  # Iyy
        I[3,3] += m * (r[1]^2 + r[2]^2)  # Izz
        
        # Products of inertia
        I[1,2] = I[2,1] -= m * r[1] * r[2]  # Ixy
        I[1,3] = I[3,1] -= m * r[1] * r[3]  # Ixz
        I[2,3] = I[3,2] -= m * r[2] * r[3]  # Iyz
    end
    
    # Find principal axes
    eigenvals, eigenvecs = eigen(I)
    
    # Sort by magnitude
    p = sortperm(eigenvals)
    eigenvals = eigenvals[p]
    eigenvecs = eigenvecs[:, p]
    
    # Ensure right-handed coordinate system
    if det(eigenvecs) < 0
        eigenvecs[:, 3] *= -1
    end

    # Store results
    s.I_kite .= eigenvals
    s.R_b_p .= eigenvecs
    s.Q_p_b .= rotation_matrix_to_quaternion(s.R_b_p')
    return nothing
end

function calc_pos_principal!(s::KPSQ)
    # pos in principal frame relative to com
    mass_per_area = s.set.mass / ((s.set.middle_length + s.set.tip_length) * 0.5 * s.set.width)
    n = s.set.aero_surfaces
    s.α_l       = π/2 - s.set.width/2/s.set.radius
    α_middle    = π/2
    dα          = (α_middle - s.α_l) / n
    for i in 1:2n
        if i <= n
            α = s.α_l + -dα/2 + i * dα
            kite_length = s.set.tip_length + (s.set.middle_length-s.set.tip_length) * (α - s.α_l) / (π/2 - s.α_l) # TODO: kite length gets less with flap turning
        else
            α = pi - (s.α_l + -dα/2 + (i-n) * dα)
            kite_length = s.set.tip_length + (s.set.middle_length-s.set.tip_length) * (π - s.α_l - α) / (π/2 - s.α_l)
        end
        area = kite_length * s.set.width/(2n)
        s.seg_mass[i] = area * mass_per_area

        s.seg_com_pos_p[:, i] .= s.R_b_p * ([-0.5 * kite_length, cos(α) * s.set.radius, sin(α) * s.set.radius] .+ s.pos_circle_center_b)
        s.seg_cop_pos_b[:, i] .= [-0.75 * kite_length, cos(α) * s.set.radius, sin(α) * s.set.radius] + s.pos_circle_center_b
        s.seg_cop_pos_p[:, i] .= s.R_b_p * s.seg_cop_pos_b[:, i]
    end
    s.pos_C_b = s.pos_circle_center_b .+ [-0.75s.kite_length_D, 0.0, s.set.radius - s.set.bridle_center_distance]
    s.pos_C_p .= s.R_b_p * s.pos_C_b
    s.pos_A_b .= [0.0, cos(s.α_D), 0.0] .+ s.pos_C_b
    s.pos_B_b .= [0.0, -cos(s.α_D), 0.0] .+ s.pos_C_b
    return nothing
end

function init_pos!(s::KPSQ; distance=s.measure.tether_length[3])
    # ground points
    s.pos .= 0.0
    s.kite_pos .= 0.0

    # init last tether points
    s.pos[:, s.i_A] .= rotate_around_z(rotate_around_y([distance, 0, 0], -s.measure.elevation_left), s.measure.azimuth_left)
    s.pos[:, s.i_B] .= rotate_around_z(rotate_around_y([distance, 0, 0], -s.measure.elevation_right), s.measure.azimuth_right)
    s.pos[:, s.i_C] .= 0.5s.pos[:, s.i_A] + 0.5s.pos[:, s.i_B]
    s.e_y .= normalize(s.pos[:, s.i_A] - s.pos[:, s.i_B])

    # init middle tether
    angular_acc = s.measure.tether_acc / s.set.drum_radius
    net_torque = angular_acc * s.set.inertia_total # TODO: check if inertia is correct
    tether_force = (net_torque - s.measure.set_values) / s.set.drum_radius
    s.pos[:, 3:3:s.i_C] .= calc_expected_pos_vel(s, s.pos[:, s.i_C][1], s.pos[:, s.i_C][2], s.pos[:, s.i_C][3], 
        0, 0, s.measure.tether_length[3], tether_force[3], s.c_spring[3])[1, :, :]
    s.e_z .= normalize(s.pos[:, s.i_C] - s.pos[:, s.i_C-3])
    s.e_x .= s.e_y × s.e_z
    R_b_w = hcat(s.e_x, s.e_y, s.e_z)
    s.Q_p_w .= rotation_matrix_to_quaternion(R_b_w * s.R_b_p')
    s.kite_pos .= s.pos[:, s.i_C] - R_b_w * s.pos_C_b
    
    # init tether connection points
    s.pos_D_b .= [0.0, cos(s.α_D) * s.set.radius, sin(s.α_D) * s.set.radius]
    s.pos_E_b .= [0.0, -cos(s.α_D) * s.set.radius, sin(s.α_D) * s.set.radius]
    e_r_D = -normalize(s.pos_D_b)
    e_r_E = -normalize(s.pos_E_b)
    te_length = s.kite_length_D/4
    angle_te_c = 0.0
    angle_te_d = 0.0
    s.pos[:, s.i_A] .= s.pos[:, s.i_C] + s.e_y * s.pos_D_b[2] + s.e_x * te_length * cos(angle_te_c) + e_r_D * te_length * sin(angle_te_c)
    s.pos[:, s.i_B] .= s.pos[:, s.i_C] + s.e_y * s.pos_E_b[2] + s.e_x * te_length * cos(angle_te_d) + e_r_E * te_length * sin(angle_te_d)

    # init left and right tether
    s.pos[:, 1:3:s.i_A] .= calc_expected_pos_vel(s, s.pos[:, s.i_A][1], s.pos[:, s.i_A][2], s.pos[:, s.i_A][3], 
        0, 0, s.measure.tether_length[1], tether_force[1], s.c_spring[1])[1, :, :]
    s.pos[:, 2:3:s.i_B] .= calc_expected_pos_vel(s, s.pos[:, s.i_B][1], s.pos[:, s.i_B][2], s.pos[:, s.i_B][3], 
        0, 0, s.measure.tether_length[2], tether_force[2], s.c_spring[2])[1, :, :]
    return s.pos
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

function init_inner(s::KPSQ; delta=0.0)
    pos_, vel_, acc_ = init_pos_vel_acc(s; delta=delta)
    # remove last left and right tether point and replace them by the length from C and D
    pos = vcat(
        pos_[4:s.i_A-1], 
        pos_[s.i_E:end],
    )
    len = vcat( # connection length
        norm(pos_[s.i_C]-pos_[s.i_A]), # left line connection distance
        norm(pos_[s.i_D]-pos_[s.i_B]), # right line connection distance
    )
    vel = vcat(
        vel_[4:s.i_A-1], 
        vel_[s.i_E:end],
    )
    acc = vcat(
        acc_[4:s.i_A-1],
        acc_[s.i_E:end],
    )
    vcat(pos, vel, len, [0,0]), vcat(vel, acc, [0,0], [0,0])
end

# same as above, but returns a tuple of two one dimensional arrays
function init(s::KPS4, X=zeros(2 * (s.set.segments+KITE_PARTICLES-1)+1); old=false, delta=0.0, upwind_dir=nothing)
    res1_, res2_ = init_inner(s, X; old=old, delta = delta)
    res1__ = reduce(vcat, res1_)
    res2__ = reduce(vcat, res2_)
    if !isnothing(upwind_dir)
        res1__ = turn(res1__, upwind_dir)
        res2__ = turn(res2__, upwind_dir)
    end
    res1, res2  = vcat(res1__, [s.l_tether, 0]),  vcat(res2__,[0,0])
    MVector{6*(s.set.segments+KITE_PARTICLES)+2, SimFloat}(res1), MVector{6*(s.set.segments+KITE_PARTICLES)+2, SimFloat}(res2)
end

function init(s::KPSQ; delta=0.0)
    y_, yd_ = init_inner(s; delta = delta)
    y = vcat(reduce(vcat, y_), reduce(vcat,[s.tether_lengths, zeros(3)]))
    yd = vcat(reduce(vcat, yd_), zeros(6))
    MVector{6*(s.i_A-5)+4+6, SimFloat}(y), MVector{6*(s.i_A-5)+4+6, SimFloat}(yd)
end

# rotate a 3d vector around the y axis in the xz plane - following the right hand rule
function rotate_around_y(vec, angle)
    result = similar(vec)
    result[1] = cos(angle) * vec[1] + sin(angle) * vec[3]
    result[2] = vec[2]
    result[3] = -sin(angle) * vec[1] + cos(angle) * vec[3]
    result
end

# rotate a 3d vector around the z axis in the yx plane - following the right hand rule
function rotate_around_z(vec, angle)
    result = similar(vec)
    result[1] = cos(angle) * vec[1] - sin(angle) * vec[2]
    result[2] = sin(angle) * vec[1] + cos(angle) * vec[2]
    result[3] = vec[3]
    result
end

function rotate_v_around_k(v, k, θ)
    k = normalize(k)
    v_rot = v * cos(θ) + (k × v) * sin(θ)  + k * (k ⋅ v) * (1 - cos(θ))
    return v_rot
end