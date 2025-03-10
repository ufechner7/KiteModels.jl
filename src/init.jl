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
    vel = zeros(MVector{s.set.segments+1+KITE_PARTICLES, KVec3})
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
    for i in 2:s.set.segments+1
        vel[i] .+= (pos[i+1] - pos[i]) * (s.set.v_reel_out*(i-1)/s.set.segments)
    end
    # the velocity vector of the kite particles is the same as the velocity of the last tether pointj
    for i in s.set.segments+2:s.set.segments+KITE_PARTICLES+1
        vel[i] .+= vel[s.set.segments+1] 
    end
    pos, vel, acc
end

function set_initial_velocity!(s::KPS4)
    for i in 2:s.set.segments+1
        s.vel[i] .= (s.pos[i+1] - s.pos[i]) * (s.set.v_reel_out*(i-1)/s.set.segments)
    end
    # the velocity vector of the kite particles is the same as the velocity of the last tether point
    for i in s.set.segments+2:s.set.segments+KITE_PARTICLES+1
        s.vel[i] .= s.vel[s.set.segments+1] 
    end
end

# function find_bridle_gammas!(s::KPSQ, wing::KiteWing; n_groups=4)
#     total_area = wing.area_interp(0.0)
#     function equations!(F, gamma, p)
#         F[1] = wing.area_interp(gamma[1]) - 1/4 * total_area
#         F[2] = wing.area_interp(gamma[2]) - 3/4 * total_area
#     end
#     gamma_tip = wing.gamma_tip
#     gamma0 = [-gamma_tip + 0.25*gamma_tip, -gamma_tip + 0.75*gamma_tip]
#     prob = NonlinearProblem(equations!, gamma0, nothing)
    
#     result = NonlinearSolve.solve(prob, NewtonRaphson())

#     bridle_gamma = zeros(n_groups)
#     bridle_gamma[1:2] .= result.u
#     bridle_gamma[3:4] .= -bridle_gamma[1:2]
#     return bridle_gamma
# end

function find_bridle_gammas!(s::KPSQ, wing::KiteWing; n_groups=4)
    @assert iseven(n_groups) "Number of groups must be even"
    half_area = wing.area_interp(0.0)
    
    # Create target areas - evenly distributed on left side of the wing
    target_areas = [(i / n_groups) * (half_area) for i in 1:n_groups-1]
    
    function equations!(F, gamma, p)
        for i in 1:n_groups-1
            F[i] = wing.area_interp(gamma[i]) - target_areas[i]
        end
    end
    
    # Initial guess: evenly spaced between -gamma_tip and middle
    gamma_tip = wing.gamma_tip
    gamma0 = [-gamma_tip + (i / n_groups) * gamma_tip for i in 1:n_groups-1]
    
    prob = NonlinearProblem(equations!, gamma0, nothing)
    result = NonlinearSolve.solve(prob, NewtonRaphson())
    
    # Mirror the solution for both sides
    bridle_gamma = zeros(n_groups)
    bridle_gamma[1:n_groups÷2] .= result.u[1:2:end]
    bridle_gamma[n_groups÷2+1:end] .= -reverse(result.u[1:2:end])

    limits = [zeros(2) for _ in 1:n_groups]
    for i in eachindex(limits[1:n_groups÷2])
        if i == 1
            limits[i] .= (-gamma_tip, result.u[2i])
        elseif 2i == length(result.u)+1
            limits[i] .= (result.u[2i-2], 0.0)
        else
            limits[i] .= (result.u[2i-2], result.u[2i])
        end
        limits[n_groups+1-i] .= -reverse(limits[i])
    end
    return bridle_gamma, limits
end

function init_bridle_pos!(s)
    bridle_pos_b = zeros(SimFloat, 3, 4, 4) # xyz, length, width
    bridle_gamma = zeros(SimFloat, 4)

    calc_inertia!(s, s.wing)
    find_attachment_gamma!(s, bridle_gamma)

    bridle_fracs = (s.set.bridle_connect[2:5] .- s.set.bridle_connect[1]) / (s.set.bridle_connect[1] - s.set.bridle_connect[5])
    @show bridle_fracs
    for (i, gamma) in enumerate(bridle_gamma)
        for (j, frac) in enumerate(bridle_fracs)
            bridle_pos_b[:, j, i] .= s.pos_circle_center_b + 
                [s.leading_edge(gamma) + frac * (s.trailing_edge(gamma) - s.leading_edge(gamma)), cos(gamma) * s.set.radius, sin(gamma) * s.set.radius]
        end
    end
    return bridle_pos_b
end

function calc_inertia!(s::KPSQ, wing::KiteWing)
    s.I_b .= [wing.inertia_tensor[1,1], wing.inertia_tensor[2,2], wing.inertia_tensor[3,3]]
    
    # Find principal axes
    eigenvals, eigenvecs = eigen(I)
    
    # Sort by magnitude
    p = sortperm(eigenvals)
    eigenvals = eigenvals[p]
    eigenvecs = eigenvecs[:, p]
    
    # Ensure right-handed coordinate system
    if det(eigenvecs) < 0
        eigenvecs[:, 3] .*= -1
    end

    # Store results
    s.I_p .= eigenvals
    s.R_b_p .= eigenvecs
    s.Q_p_b .= rotation_matrix_to_quaternion(s.R_b_p')
    return nothing
end

function calc_pos_principal!(s::KPSQ)
    # pos in principal frame relative to com
    mass_per_area = s.set.mass / ((s.set.middle_length + s.set.tip_length) * 0.5 * s.set.width)
    n = s.set.aero_surfaces
    s.gamma_l       = π/2 - s.set.width/2/s.set.radius
    gamma_middle    = π/2
    dgamma          = (gamma_middle - s.gamma_l) / n
    for i in 1:2n
        if i <= n
            gamma = s.gamma_l + -dgamma/2 + i * dgamma
            kite_length = s.set.tip_length + (s.set.middle_length-s.set.tip_length) * (gamma - s.gamma_l) / (π/2 - s.gamma_l)
        else
            gamma = pi - (s.gamma_l + -dgamma/2 + (i-n) * dgamma)
            kite_length = s.set.tip_length + (s.set.middle_length-s.set.tip_length) * (π - s.gamma_l - gamma) / (π/2 - s.gamma_l)
        end
        area = kite_length * s.set.width/(2n)
        s.seg_mass[i] = area * mass_per_area

        s.seg_com_pos_p[:, i] .= s.R_b_p * ([-0.5 * kite_length, cos(gamma) * s.set.radius, sin(gamma) * s.set.radius] .+ s.pos_circle_center_b)
        s.seg_cop_pos_b[:, i] .= [-0.75 * kite_length, cos(gamma) * s.set.radius, sin(gamma) * s.set.radius] + s.pos_circle_center_b
        s.seg_cop_pos_p[:, i] .= s.R_b_p * s.seg_cop_pos_b[:, i]
    end
    s.pos_C_b = s.pos_circle_center_b .+ [-0.75s.kite_length_D, 0.0, s.set.radius - s.set.bridle_center_distance]
    s.pos_C_p .= s.R_b_p * s.pos_C_b
    s.pos_A_b .= [0.0, cos(s.gamma_D), 0.0] .+ s.pos_C_b
    s.pos_B_b .= [0.0, -cos(s.gamma_D), 0.0] .+ s.pos_C_b
    return nothing
end

function init_pos!(s::KPSQ; distance = nothing, te_angle = s.te_angle, kite_angle = 0.0)
    # ground points
    s.pos .= 0.0
    s.kite_pos .= 0.0

    # init last tether points
    angular_acc = s.measure.tether_acc / s.set.drum_radius
    net_torque = angular_acc * s.set.inertia_total
    tether_force = (net_torque - s.measure.set_values) / s.set.drum_radius
    if isnothing(distance)
        distance = s.measure.tether_length[3] + tether_force[3] / (s.c_spring[3]/s.measure.tether_length[3]) - 1e-3
    end
    s.pos[:, s.i_A] .= rotate_around_z(rotate_around_y([distance, 0, 0], -s.measure.sphere_pos[1, 1]), s.measure.sphere_pos[2, 1])
    s.pos[:, s.i_B] .= rotate_around_z(rotate_around_y([distance, 0, 0], -s.measure.sphere_pos[1, 2]), s.measure.sphere_pos[2, 2])
    s.pos[:, s.i_C] .= 0.5 .* s.pos[:, s.i_A] .+ 0.5 .* s.pos[:, s.i_B]
    s.e_y .= normalize(s.pos[:, s.i_A] .- s.pos[:, s.i_B])

    # init middle tether
    s.pos[:, 3:3:s.i_C] .= calc_expected_pos_vel(s, s.pos[:, s.i_C], 
        0, 0, s.measure.tether_length[3], tether_force[3], s.c_spring[3])[1, :, :]
    s.e_z .= normalize(s.pos[:, s.i_C] - s.pos[:, s.i_C-3])
    s.e_z .= rotate_v_around_k(s.e_z, s.e_y, kite_angle)
    s.e_x .= s.e_y × s.e_z
    R_b_w = hcat(s.e_x, s.e_y, s.e_z)
    s.Q_p_w .= rotation_matrix_to_quaternion(R_b_w * s.R_b_p')
    s.kite_pos .= s.pos[:, s.i_C] - R_b_w * s.pos_C_b
    
    # init tether connection points
    s.pos_D_b .= [0.0, cos(s.gamma_D) * s.set.radius, sin(s.gamma_D) * s.set.radius]
    s.pos_E_b .= [0.0, -cos(s.gamma_D) * s.set.radius, sin(s.gamma_D) * s.set.radius]
    e_r_D = -normalize(s.pos_D_b)
    e_r_E = -normalize(s.pos_E_b)
    te_length = s.kite_length_D/4
    s.pos[:, s.i_A] .= s.pos[:, s.i_C] .+ s.e_y .* s.pos_D_b[2] .+ s.e_x .* te_length * cos(te_angle[1]) .+ e_r_D .* te_length * sin(te_angle[1])
    s.pos[:, s.i_B] .= s.pos[:, s.i_C] .+ s.e_y .* s.pos_E_b[2] .+ s.e_x .* te_length * cos(te_angle[2]) .+ e_r_E .* te_length * sin(te_angle[2])

    # init left and right tether
    s.pos[:, 1:3:s.i_A] .= calc_expected_pos_vel(s, s.pos[:, s.i_A], 
        0, 0, s.measure.tether_length[1], tether_force[1], s.c_spring[1])[1, :, :]
    s.pos[:, 2:3:s.i_B] .= calc_expected_pos_vel(s, s.pos[:, s.i_B], 
        0, 0, s.measure.tether_length[2], tether_force[2], s.c_spring[2])[1, :, :]
    return s.pos
end

function reinit!(s::KPSQ; init=false, new_sys=true, initial_damping=0., damping_time=10)
    dt = 40.
    tspan = (0., dt)
    if new_sys
        solver = QBDF( # https://docs.sciml.ai/SciMLBenchmarksOutput/stable/#Results
            autodiff=AutoFiniteDiff()
        )
        if !init
            isys, u0map, p0map = model!(s; init=true)
            s.init_prob = ODEProblem(isys, u0map, tspan, p0map)
            isys = s.init_prob.f.sys
            s.init_integrator = OrdinaryDiffEqCore.init(s.init_prob, solver; dt, abstol=s.set.abs_tol, reltol=s.set.rel_tol, save_on=false)
        else
            isys = s.prob.f.sys
            s.init_prob = s.prob
            s.init_integrator = s.integrator
        end

        set_initial_measure = setp(s.init_prob, (
            isys.measured_wind_dir_gnd,
            isys.measured_sphere_pos,
            isys.measured_sphere_vel,
            isys.measured_sphere_acc,
            isys.measured_tether_length,
            isys.measured_tether_vel,
            isys.measured_tether_acc,
        ))
        set_initial_set_values = setu(isys, isys.set_values)
        
        s.set_initial_measure = (prob, val) -> set_initial_measure(prob, val)
        s.set_initial_set_values = (prob, val) -> set_initial_set_values(prob, val)
        
        diff_state = unknowns(s.prob.f.sys)
        get_diff_state = getu(isys, diff_state)
        set_diff_state = setu(s.prob, diff_state)
        s.get_diff_state = (prob) -> get_diff_state(prob)
        s.set_diff_state = (prob, val) -> set_diff_state(prob, val)
    end
    s.set_initial_measure(s.init_prob, (
        s.measure.wind_dir_gnd,
        s.measure.sphere_pos,
        s.measure.sphere_vel,
        s.measure.sphere_acc,
        s.measure.tether_length,
        s.measure.tether_vel,
        s.measure.tether_acc,
    ))
    s.set_initial_set_values(s.init_prob, s.measure.set_values)
    s.init_prob = remake(s.init_prob)
    OrdinaryDiffEqCore.reinit!(s.init_integrator, s.init_prob.u0)
    OrdinaryDiffEqCore.step!(s.init_integrator, dt, true) # TODO: step until a certain amount of time has passed
    diff_state = s.get_diff_state(s.init_integrator)
    s.set_diff_state(s.prob, diff_state)
    s.prob = remake(s.prob)
    OrdinaryDiffEqCore.reinit!(s.integrator, s.prob.u0)
    nothing
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
    set_initial_velocity!(s)
    for i in 1:Int(length(res1__)/6)
        j = i + s.set.segments+KITE_PARTICLES
        # println("i, j:", i, ", ", j)
        res1__[3*(j-1)+1:3*j] .= s.vel[i+1]
        res2__[3*(i-1)+1:3*i] .= s.vel[i+1]
    end
    res1, res2  = vcat(res1__, [s.l_tether, s.set.v_reel_out]),  vcat(res2__,[0,0])
    MVector{6*(s.set.segments+KITE_PARTICLES)+2, SimFloat}(res1), MVector{6*(s.set.segments+KITE_PARTICLES)+2, SimFloat}(res2)
end

function init(s::KPSQ; delta=0.0)
    y_, yd_ = init_inner(s; delta = delta)
    y = vcat(reduce(vcat, y_), reduce(vcat,[s.tether_lengths, zeros(3)]))
    yd = vcat(reduce(vcat, yd_), zeros(6))
    MVector{6*(s.i_A-5)+4+6, SimFloat}(y), MVector{6*(s.i_A-5)+4+6, SimFloat}(yd)
end

# rotate a 3d vector around the y axis in the xz plane - following the right hand rule
function rotate_around_y(vec, angle::T) where T
    result = zeros(T, 3)
    result[1] = cos(angle) * vec[1] + sin(angle) * vec[3]
    result[2] = vec[2]
    result[3] = -sin(angle) * vec[1] + cos(angle) * vec[3]
    result
end

# rotate a 3d vector around the z axis in the yx plane - following the right hand rule
function rotate_around_z(vec, angle::T) where T
    result = zeros(T, 3)
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