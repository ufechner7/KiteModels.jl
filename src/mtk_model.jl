# ==================== mtk model functions ================================================
# Implementation of the three-line model using ModellingToolkit.jl

function smooth_sign_ϵ(x; ϵ = 1e-3)
    return x / √(x^2 + ϵ^2)
end
@register_symbolic smooth_sign_ϵ(x)

function calc_acc_speed(motor::AsyncMachine, tether_vel, norm_, set_speed)
    calc_acceleration(motor, tether_vel, norm_; set_speed, set_torque=nothing, use_brake=false) # TODO: add brake setting
end
@register_symbolic calc_acc_speed(motor::AsyncMachine, tether_vel, norm_, set_speed)

function calc_acc_torque(motor::TorqueControlledMachine, tether_vel, norm_, set_torque)
    calc_acceleration(motor, tether_vel, norm_; set_speed=nothing, set_torque, use_brake=false)
end
@register_symbolic calc_acc_torque(motor::TorqueControlledMachine, tether_vel, norm_, set_torque)

function sym_interp(interp::Function, aoa, trailing_edge_angle)
    return interp(rad2deg(aoa), rad2deg(trailing_edge_angle-aoa)) # TODO: register callable struct https://docs.sciml.ai/Symbolics/dev/manual/functions/#Symbolics.@register_array_symbolic
end
@register_symbolic sym_interp(interp::Function, aoa, trailing_edge_angle)

function normalize(vec)
    return vec / norm(vec)
end
@register_symbolic normalize(vec)

function convert_pos_vel(s::KPSQ, pos_, vel_)
    pos = Array{Union{Nothing, Float64}}(nothing, 3, s.i_A)
    vel = Array{Union{Nothing, Float64}}(nothing, 3, s.i_A)
    [pos[:,i] .= pos_[i] for i in 1:s.i_A-1]
    [vel[:,i] .= vel_[i] for i in 1:s.i_A-1]
    [pos[:,i] .= pos_[i] for i in s.i_B+1:s.i_A]
    [vel[:,i] .= vel_[i] for i in s.i_B+1:s.i_A]
    return pos, vel
end

function rotate_by_quaternion(q, v)
    p = [0, v...]
    return quaternion_multiply(q, quaternion_multiply(p, quaternion_conjugate(q)))[2:4]
end

function quaternion_conjugate(q)
    return [q[1], -q[2], -q[3], -q[4]]
end

function quaternion_multiply(q1, q2)
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    w = w1*w2 - x1*x2 - y1*y2 - z1*z2
    x = w1*x2 + x1*w2 + y1*z2 - z1*y2
    y = w1*y2 - x1*z2 + y1*w2 + z1*x2
    z = w1*z2 + x1*y2 - y1*x2 + z1*w2
    return [w, x, y, z]
end

function quaternion_to_rotation_matrix(q)
    w, x, y, z = q[1], q[2], q[3], q[4]
    
    return [
        1 - 2*(y*y + z*z)  2*(x*y - z*w)      2*(x*z + y*w);
        2*(x*y + z*w)      1 - 2*(x*x + z*z)  2*(y*z - x*w);
        2*(x*z - y*w)      2*(y*z + x*w)      1 - 2*(x*x + y*y)
    ]
end

function rotation_matrix_to_quaternion(R)
    tr_ = tr(R)
    
    if tr_ > 0
        S = sqrt(tr_ + 1.0) * 2
        w = 0.25 * S
        x = (R[3,2] - R[2,3]) / S
        y = (R[1,3] - R[3,1]) / S
        z = (R[2,1] - R[1,2]) / S
    elseif (R[1,1] > R[2,2]) && (R[1,1] > R[3,3])
        S = sqrt(1.0 + R[1,1] - R[2,2] - R[3,3]) * 2
        w = (R[3,2] - R[2,3]) / S
        x = 0.25 * S
        y = (R[1,2] + R[2,1]) / S
        z = (R[1,3] + R[3,1]) / S
    elseif R[2,2] > R[3,3]
        S = sqrt(1.0 + R[2,2] - R[1,1] - R[3,3]) * 2
        w = (R[1,3] - R[3,1]) / S
        x = (R[1,2] + R[2,1]) / S
        y = 0.25 * S
        z = (R[2,3] + R[3,2]) / S
    else
        S = sqrt(1.0 + R[3,3] - R[1,1] - R[2,2]) * 2
        w = (R[2,1] - R[1,2]) / S
        x = (R[1,3] + R[3,1]) / S
        y = (R[2,3] + R[3,2]) / S
        z = 0.25 * S
    end
    
    return [w, x, y, z]
end

"""
Calculate the forces acting on the kite inertia particle.
"""
function calc_kite_forces!(s::KPSQ, seqs, force_eqs, force, torque_p, R_b_w, kite_vel, kite_acc, ω_b, t, e_x, e_z, rho, v_wind, trailing_edge_angle)
    n = s.set.aero_surfaces

    # integrating loop variables, iterating over 2n segments
    @variables begin
        aero_pos(t)[1:3, 1:2n]
        e_r(t)[1:3, 1:2n]
        e_te(t)[1:3, 1:2n] # clockwise trailing edge of flap vector
        y_l(t)[1:2n]
        v_a(t)[1:3, 1:2n]
        e_drift(t)[1:3, 1:2n]
        v_a_xr(t)[1:3, 1:2n]
        aoa(t)[1:n*2]
        seg_cl(t)[1:n*2]
        seg_cd(t)[1:n*2]
        seg_L(t)[1:3, 1:2n]
        seg_D(t)[1:3, 1:2n]
        seg_aero_force(t)[1:3, 1:2n]
        total_kite_force(t)[1:3]
        seg_g_force(t)[1:3, 1:2n]
        seg_te_force(t)[1:3, 1:2n]
        seg_trailing_edge_angle(t)[1:2n] # flap angle relative to -e_x, e_z
        ram_force(t)[1:2n]
        te_force(t)[1:2n]
        seg_aero_torque_p(t)[1:3, 1:2n]
        seg_gravity_torque_p(t)[1:3, 1:2n]
        C_torque_p(t)[1:3] # torque on point C in principal reference frame
        seg_vel(t)[1:3, 1:2n]
        F(t)[1:3, 1:2n]
        r(t)[1:3, 1:2n]
        L_C(t)[1:3]
        L_D(t)[1:3]
        D_C(t)[1:3]
        D_D(t)[1:3]
        F_te_C(t)[1:3] # trailing edge clockwise force
        F_te_D(t)[1:3]
    end

    s.α_l         = π/2 - s.set.width/2/s.set.radius
    α_middle    = π/2
    dα          = (α_middle - s.α_l) / n
    α_m_l = π/2 - s.set.min_steering_line_distance/s.set.radius # end of steering lines on left side
    α_m_r = π/2 + s.set.min_steering_line_distance/s.set.radius # end of steering lines on right side
    
    for i in 1:2n
        if i <= n
            α = s.α_l + -dα/2 + i * dα
            kite_length = s.set.tip_length + (s.set.middle_length-s.set.tip_length) * (α - s.α_l) / (π/2 - s.α_l) # TODO: kite length gets less with flap turning
        else
            α = pi - (s.α_l + -dα/2 + (i-n) * dα)
            kite_length = s.set.tip_length + (s.set.middle_length-s.set.tip_length) * (π - s.α_l - α) / (π/2 - s.α_l)
        end
        seg_flap_height = kite_length * s.set.flap_height

        seqs = [
            seqs
            e_r[:, i]       ~ rotate_v_around_k(e_z, e_x, 0.5π + α)
            seg_vel[:, i]   ~ kite_vel + R_b_w * (ω_b × s.seg_cop_pos_b[:, i])
            v_a[:, i]       ~ v_wind .- seg_vel[:, i]
            e_drift[:, i]   ~ (e_x × e_r[:, i])
            v_a_xr[:, i]    ~ v_a[:, i] .- (v_a[:, i] ⋅ e_drift[:, i]) .* e_drift[:, i]

            α < α_m_l ?
                seg_trailing_edge_angle[i] ~ trailing_edge_angle[1]  :
            α > α_m_r ?
                seg_trailing_edge_angle[i] ~ trailing_edge_angle[2] :
                seg_trailing_edge_angle[i] ~ ((trailing_edge_angle[2] - trailing_edge_angle[1]) / (α_m_r - α_m_l) * (α - α_m_l) + trailing_edge_angle[1])

            aoa[i]      ~ -asin2(normalize(v_a_xr[:, i]) ⋅ e_r[:, i]) + deg2rad(s.set.alpha_zero)
            seg_cl[i]   ~ sym_interp(s.cl_interp, aoa[i], seg_trailing_edge_angle[i])
            seg_cd[i]   ~ sym_interp(s.cd_interp, aoa[i], seg_trailing_edge_angle[i])

            seg_L[:, i] ~ 0.5 * rho * (norm(v_a_xr[:, i]))^2 * s.set.radius * dα * kite_length * seg_cl[i] * 
                                normalize(v_a_xr[:, i] × e_drift[:, i])
            seg_D[:, i] ~ 0.5 * rho * norm(v_a_xr[:, i]) * s.set.radius * dα * kite_length * seg_cd[i] *
                                v_a_xr[:, i]

            e_te[:, i] ~ -e_x * sin(seg_trailing_edge_angle[i]) + e_r[:, i] * cos(seg_trailing_edge_angle[i])
            ram_force[i] ~ smooth_sign_ϵ(deg2rad(s.set.alpha_zero) - seg_trailing_edge_angle[i]; s.ϵ) *
                        rho * norm(v_a[:, i])^2 * seg_flap_height * s.set.radius * dα * (seg_flap_height/2) / (kite_length/4)
            te_force[i] ~ 0.5 * rho * (norm(v_a_xr[:, i]))^2 * s.set.radius * dα * kite_length * 
                        sym_interp(s.c_te_interp, aoa[i], seg_trailing_edge_angle[i])
            seg_te_force[:, i] ~ (ram_force[i] + te_force[i]) * e_te[:, i]
            
            seg_aero_force[:, i] ~ seg_L[:, i] + seg_D[:, i]
            seg_g_force[:, i] ~ [0.0, 0.0, -G_EARTH * s.seg_mass[i]]
            seg_aero_torque_p[:, i] ~ s.seg_cop_pos_p[:, i] × (s.R_b_p * R_b_w' * seg_aero_force[:, i])
            seg_gravity_torque_p[:, i] ~ s.seg_com_pos_p[:, i] × (s.R_b_p * R_b_w' * seg_g_force[:, i])
        ]
    end

    seqs = [
        seqs
        total_kite_force ~ [sum(seg_aero_force[i, :]) + sum(seg_g_force[i, :]) + force[i, s.i_C] for i in 1:3]
        kite_acc ~ total_kite_force / s.set.mass
        C_torque_p ~ s.pos_C_p × (s.R_b_p * R_b_w' * force[:, s.i_C])
        torque_p ~ [sum(seg_aero_torque_p[i, :]) + sum(seg_gravity_torque_p[i, :]) + C_torque_p[i] for i in 1:3]
        F_te_C ~ [sum(seg_te_force[i, 1:n]) for i in 1:3]
        F_te_D ~ [sum(seg_te_force[i, n+1:2n]) for i in 1:3]
    ]

    # longtitudinal force
    # F_inside_flap = P * A
    # F_inside_flap = rho * norm(v)^2 * flap_height * width
    # F_inside_flap = rho * norm(v)^2 * flap_height * radius * dα
    # F_trailing_edge = -F_inside_flap * (flap_height/2) / (kite_length/4) if trailing_edge_angle > 0 clockwise force
    # F_trailing_edge = F_inside_flap * (flap_height/2) / (kite_length/4) if trailing_edge_angle < 0 clockwise force
    # dF_te_dα = rho * norm(v)^2 * flap_height * radius
    # flap_height = height_middle * kite_length / middle_length
    
    force_eqs[:,s.i_A] .= (force[:, s.i_A] .~ F_te_C + [0.0, 0.0, -G_EARTH*(s.set.mass/8)])
    force_eqs[:,s.i_B] .= (force[:, s.i_B] .~ F_te_D + [0.0, 0.0, -G_EARTH*(s.set.mass/8)])
    return seqs, force_eqs
end

""" 
    calc_particle_forces!(s::KPSQ, seqs, force_eqs, force, pos1, pos2, vel1, vel2, length, c_spring, damping, 
                          rho, i, l_0, k, c, segment, rel_vel, av_vel, norm1, unit_vector, k2, c1, spring_vel,
                          spring_force, v_apparent, v_wind_tether, area, v_app_perp, half_drag_force)

Calculate the drag force and spring force of the tether segment, defined by the parameters pos1, pos2, vel1 and vel2
and distribute it equally on the two particles, that are attached to the segment.
The result is stored in the array s.forces. 
"""
function calc_particle_forces!(s::KPSQ, seqs, force_eqs, force, p1, p2, pos1, pos2, vel1, vel2, length, c_spring, 
    damping, rho, i, l_0, k, c, segment, rel_vel, av_vel, norm1, unit_vector, k2, c1, c2, spring_vel, perp_vel,
            spring_force, v_apparent, v_wind_tether, area, v_app_perp, half_drag_force)
    d_tether = s.set.d_tether/1000.0
    seqs = [
        seqs
        l_0 ~ length[(i-1) % 3 + 1] # Unstressed length
        k   ~ c_spring[(i-1) % 3 + 1] # Spring constant
        c   ~ damping[(i-1) % 3 + 1] # Damping coefficient    
        segment     .~ pos1 - pos2
        rel_vel     .~ vel1 - vel2
        av_vel      .~ 0.5 * (vel1 + vel2)
        norm1        ~ norm(segment)
        unit_vector .~ segment / norm1
        k2           ~ 0.1 * k  # compression stiffness tether segments
        c1           ~ 6.0 * c  # damping kite segments
        c2           ~ 0.05 * c  # damping perpendicular
        spring_vel   ~ rel_vel ⋅ unit_vector
        perp_vel    .~ rel_vel .- spring_vel * unit_vector
    ]

    for j in 1:3
        seqs = [
            seqs
            spring_force[j] ~
                ((k  * (l_0 - norm1) - c * spring_vel) * unit_vector[j]) * (1 + smooth_sign_ϵ(norm1 - l_0; s.ϵ)) / 2 +
                ((k2 * (l_0 - norm1) - c * spring_vel) * unit_vector[j]) * (1 - smooth_sign_ϵ(norm1 - l_0; s.ϵ)) / 2 -
                c2 * perp_vel[j]
        ]
    end
    seqs = [
        seqs
        v_apparent       ~ v_wind_tether - av_vel
        i >= s.i_A ?
            area             ~ norm1 * d_tether * 10 : # 10 is the number of parallel lines in the bridle system
            area             ~ norm1 * d_tether * (1 + (i%3 == 0)) # double area for middle tether
        v_app_perp       ~ v_apparent - (v_apparent ⋅ unit_vector) * unit_vector
        half_drag_force .~ (0.25 * rho * s.set.cd_tether * norm(v_app_perp) * area) .* v_app_perp
    ]

    for j in 1:3
        force_eqs[j, p1] = 
            (force[j, p1] ~ force_eqs[j, p1].rhs + (half_drag_force[j] + spring_force[j]))
        force_eqs[j, p2] = 
            (force[j, p2] ~ force_eqs[j, p2].rhs + (half_drag_force[j] - spring_force[j]))
    end
    
    return seqs, force_eqs
end



"""
Calculate the forces, acting on all tether particles.
"""
@inline function calc_tether_forces!(s::KPSQ, seqs, force_eqs, t, force, pos, vel, length, c_spring, damping, v_wind_gnd, norm1)
    @variables begin
        height(t)[1:s.i_C-3]
        rho(t)[1:s.i_C-3]
        v_wind_tether(t)[1:3, 1:s.i_C-3]

        l_0(t)[1:s.i_C-3]
        k(t)[1:s.i_C-3]
        c(t)[1:s.i_C-3]
        segment(t)[1:3, 1:s.i_C-3]
        rel_vel(t)[1:3, 1:s.i_C-3]
        av_vel(t)[1:3, 1:s.i_C-3] 
        unit_vector(t)[1:3, 1:s.i_C-3]
        k2(t)[1:s.i_C-3]
        c1(t)[1:s.i_C-3]
        c2(t)[1:s.i_C-3]
        spring_vel(t)[1:s.i_C-3]
        perp_vel(t)[1:3, 1:s.i_C-3]
        spring_force(t)[1:3, 1:s.i_C-3]
        v_apparent(t)[1:3, 1:s.i_C-3]
        area(t)[1:s.i_C-3]
        v_app_perp(t)[1:3, 1:s.i_C-3]
        half_drag_force(t)[1:3, 1:s.i_C-3]
        gust_factor(t)
    end
    
    seqs = [seqs; D(gust_factor) ~ 0]
    for i in 1:s.i_C-3
        p1 = i  # First point nr.
        p2 = i+3
        seqs = [
            seqs
            height[i]           ~ max(0.0, 0.5 * (pos[:, p1][3] + pos[:, p2][3]))
            rho[i]              ~ calc_rho(s.am, height[i])
            v_wind_tether[:, i] ~ AtmosphericModels.calc_wind_factor(s.am, height[i], s.set.profile_law) * v_wind_gnd * abs(gust_factor)
        ]

        seqs, force_eqs = calc_particle_forces!(s, seqs, force_eqs, force, p1, p2, pos[:, p1], pos[:, p2], vel[:, p1], 
                          vel[:, p2], length, c_spring, damping, rho[i], i, l_0[i], k[i], c[i], segment[:, i], 
                          rel_vel[:, i], av_vel[:, i], norm1[i], unit_vector[:, i], k2[i], c1[i], c2[i], spring_vel[i], perp_vel[:, i],
                          spring_force[:, i], v_apparent[:, i], v_wind_tether[:, i], area[i], v_app_perp[:, i], 
                          half_drag_force[:, i])
    end

    return seqs, force_eqs
end

function expected_pos_vel(s, seqs, pos, kite_pos, kite_vel, tether_vel, tether_length, tether_force, norm1)
    @variables begin
        expected_pos(t)[1:3, 1:s.i_C]
        tether_move_vel(t)[1:3, 1:s.i_C]
        tether_kite_vel(t)[1:3, 1:s.i_C]
        expected_vel(t)[1:3, 1:s.i_C]
    end
    seqs = [
        seqs
        [expected_pos[j, 3(k-1)+i] ~ 
            calc_expected_pos_vel(s, pos[1, i+s.i_A-1], pos[2, i+s.i_A-1], pos[3, i+s.i_A-1], 
                kite_vel[i], tether_vel[i], tether_length[i], tether_force[i], s.c_spring[i])[1, j, k]
                    for i in 1:3 for j in 1:3 for k in 1:(s.i_C ÷ 3)]
        [tether_move_vel[j, 3(k-1)+i] ~ 
            calc_expected_pos_vel(s, pos[1, i+s.i_A-1], pos[2, i+s.i_A-1], pos[3, i+s.i_A-1], # TODO: A and B depend on trailing_edge_angle
                kite_vel[i], tether_vel[i], tether_length[i], tether_force[i], s.c_spring[i])[2, j, k]
                    for i in 1:3 for j in 1:3 for k in 1:(s.i_C ÷ 3)]
        [tether_kite_vel[j, i] ~ (norm(pos[:, i]) / norm(kite_pos) * 
            (kite_vel .- kite_vel ⋅ normalize(kite_pos) * normalize(kite_pos)))[j] 
                for j in 1:3 for i in 1:s.i_C]
        vec(expected_vel)                .~ vec(tether_move_vel) .+ vec(tether_kite_vel)
    ]
    return seqs
end

function scalar_eqs!(s, seqs, pos, vel, R_b_w, ω_p, ω_b, kite_pos, kite_vel, trailing_edge_angle, trailing_edge_ω, segment_length, mass_tether_particle, damping, c_spring, 
        e_y, e_z, e_x, e_r_D, e_r_E, e_te_A, e_te_B, rho_kite, tether_length, tether_vel,
        mass_per_meter, force, set_values)
    
    te_length = s.kite_length_D/4
    # last tether point vel without rotational trailing edge vel
    vel_kite_A = kite_vel + R_b_w * (ω_b × s.pos_A_b)
    vel_kite_B = kite_vel + R_b_w * (ω_b × s.pos_B_b)
    vel_kite_C = kite_vel + R_b_w * (ω_b × s.pos_C_b)
    seqs = [
        seqs
        ω_b ~ s.R_b_p' * ω_p
        
        # last tether point pos
        pos[:, s.i_C]    ~ kite_pos + R_b_w * s.pos_C_b
        pos[:, s.i_A]    ~ pos[:, s.i_C] + e_y * s.pos_D_b[2] + e_x * te_length * cos(trailing_edge_angle[1]) + # TODO: fix pos
            e_r_D * te_length * sin(trailing_edge_angle[1])
        pos[:, s.i_B]    ~ pos[:, s.i_C] + e_y * s.pos_E_b[2] + e_x * te_length * cos(trailing_edge_angle[2]) + 
            e_r_E * te_length * sin(trailing_edge_angle[2])

        # last tether point vel with rotational trailing edge vel
        vel[:, s.i_C]   ~ vel_kite_C
        vel[:, s.i_A]   ~ vel_kite_A + e_x * te_length * cos(trailing_edge_ω[1]) + e_r_D * te_length * sin(trailing_edge_ω[1])
        vel[:, s.i_B]   ~ vel_kite_B + e_x * te_length * cos(trailing_edge_ω[2]) + e_r_E * te_length * sin(trailing_edge_ω[2])

        segment_length          ~ tether_length  ./ s.set.segments
        mass_tether_particle    ~ mass_per_meter .* segment_length
        # damping                 ~ [s.damping / segment_length[1], s.damping / segment_length[2], s.damping*2 / segment_length[3]]
        damping                 ~ s.damping ./ segment_length
        # c_spring                ~ [s.c_spring / segment_length[1], s.c_spring / segment_length[2], s.c_spring*2 / segment_length[3]]
        c_spring                ~ s.c_spring ./ segment_length
        e_x     ~ R_b_w * [1, 0, 0]
        e_y     ~ R_b_w * [0, 1, 0]
        e_z     ~ R_b_w * [0, 0, 1]
        e_r_D   ~ rotate_v_around_k(e_z, e_x, 0.5π + s.α_D)
        e_r_E   ~ rotate_v_around_k(e_z, e_x, 1.5π - s.α_D)
        e_te_A  ~ -e_x * sin(trailing_edge_angle[1]) + e_r_D * cos(trailing_edge_angle[1])
        e_te_B  ~ -e_x * sin(trailing_edge_angle[2]) + e_r_E * cos(trailing_edge_angle[2])
        rho_kite        ~ calc_rho(s.am, pos[3,s.i_A])
    ]

    @variables begin
        tether_force(t)[1:3] # normalized tether forces at the winch
        heading_y(t)
        power_angle(t) # average flap angle
        power_vel(t)
        steering_angle(t) # difference between left and right flap angle
        steering_vel(t)
        tether_diff(t)
        tether_diff_vel(t)
        set_diff(t)
        azimuth(t)
        elevation(t)
        azimuth_vel(t)
        elevation_vel(t)
        distance(t)[1:s.i_C]
        distance_vel(t)
        distance_acc(t)
        kite_acc(t)[1:3]
        x_acc(t)
        y_acc(t)
        left_diff(t)
        right_diff(t)
    end
    # seqs = expected_pos_vel(s, seqs, pos, kite_pos, kite_vel, tether_vel, tether_length, tether_force, norm1)
    seqs = [
        seqs
        tether_force     ~ [norm(force[:, i]) for i in 1:3]
        heading_y       ~ calc_heading_y(-e_x)
        power_angle         ~ (trailing_edge_angle[1] + trailing_edge_angle[2]) / 2
        power_vel           ~ (trailing_edge_ω[1] + trailing_edge_ω[2]) / 2
        steering_angle      ~ trailing_edge_angle[2] - trailing_edge_angle[1]
        steering_vel        ~ trailing_edge_ω[2] - trailing_edge_ω[1]
        tether_diff         ~ tether_length[2] - tether_length[1]
        tether_diff_vel     ~ tether_vel[2] - tether_vel[1]
        set_diff            ~ set_values[2] - set_values[1]

        distance_vel        ~ kite_vel ⋅ normalize(kite_pos - pos[:, 3])
        distance_acc        ~ kite_acc ⋅ normalize(kite_pos - pos[:, 3])
        elevation           ~ atan(kite_pos[3] / kite_pos[1])
        # elevation_vel = d/dt(atan(z/x)) = (x*ż' - z*ẋ')/(x^2 + z^2) according to wolframalpha
        elevation_vel       ~ (kite_pos[1] * kite_vel[3] - kite_pos[3] * kite_vel[1]) / 
                                (kite_pos[1]^2 + kite_pos[3]^2)
        azimuth             ~ -atan(kite_pos[2] / kite_pos[1])
        # azimuth_vel = d/dt(-atan(y/x)) = (y*x' - x*y')/(x^2 + y^2)
        azimuth_vel         ~ (kite_pos[2] * kite_vel[1] - kite_pos[1] * kite_vel[2]) / 
                                (kite_pos[1]^2 + kite_pos[2]^2)
        x_acc               ~ kite_acc ⋅ e_x
        y_acc               ~ kite_acc ⋅ e_y
        left_diff           ~ tether_length[1] - tether_length[3]
        right_diff          ~ tether_length[2] - tether_length[3]
        [distance[i]        ~ norm(pos[:, i]) for i in 1:s.i_C]
    ]
    return seqs
end

function create_sys!(s::KPSQ)
    if s.torque_control
        [s.motors[i] = TorqueControlledMachine(s.set) for i in 1:3]
    else
        [s.motors[i] = AsyncMachine(s.set) for i in 1:3]
    end
    @parameters begin
        v_wind_gnd[1:3] = s.v_wind_gnd
        v_wind[1:3] = s.v_wind
    end
    @variables begin
        set_values(t)[1:3] # left right middle
        pos(t)[1:3, 1:s.i_C] # xyz pos of left right middle tether
        vel(t)[1:3, 1:s.i_C] 
        acc(t)[1:3, 1:s.i_C]
        kite_pos(t)[1:3]    # xyz position of kite in world frame
        kite_vel(t)[1:3]
        kite_acc(t)[1:3]
        R_b_w(t)[1:3, 1:3] # rotation of the kite relative to the world frame
        Q_p_w(t)[1:4] # quaternion orientation of the kite principal frame relative to the world frame
        Q_vel(t)[1:4] # quaternion rate of change
        ω_p(t)[1:3] # turn rate in principal frame
        ω_b(t)[1:3] # turn rate in body frame
        α_p(t)[1:3] # angular acceleration in principal frame
        torque_p(t)[1:3] # torque in principal frame
        trailing_edge_angle(t)[1:2]  # angle left right / C D
        trailing_edge_ω(t)[1:2]    # angular vel
        trailing_edge_α(t)[1:2]                # angular acc
        tether_length(t)[1:3]
        tether_vel(t)[1:3]
        tether_acc(t)[1:3]
        segment_length(t)[1:3]
        mass_tether_particle(t)[1:3]
        damping(t)[1:3]
        c_spring(t)[1:3]
        e_x(t)[1:3]
        e_y(t)[1:3]
        e_z(t)[1:3]
        e_r_D(t)[1:3]
        e_r_E(t)[1:3]
        e_te_A(t)[1:3]
        e_te_B(t)[1:3]
        force(t)[1:3, 1:s.i_C]
        rho_kite(t)
        norm1(t)[1:s.i_C-3]
    end
    # Collect the arrays into variables
    pos = collect(pos)
    vel = collect(vel)
    acc = collect(acc)

    deqs = []
    mass_per_meter = s.set.rho_tether * π * (s.set.d_tether/2000.0)^2
    te_length = s.kite_length_D/4 # trailing edge length
    @assert te_length != 0

    Ω = [0       -ω_p[1]   -ω_p[2]   -ω_p[3];
         ω_p[1]    0       ω_p[3]    -ω_p[2];
         ω_p[2]    -ω_p[3]   0       ω_p[1];
         ω_p[3]    ω_p[2]    -ω_p[1]   0]

    s.orient_damping = 2.0 * sqrt(maximum(s.I_kite))
    Q_b_p = quaternion_conjugate(s.Q_p_b)
    deqs = [
        deqs
        [D(Q_p_w[i]) ~ Q_vel[i] for i in 1:4]
        [Q_vel[i] ~ 0.5 * sum(Ω[i, j] * Q_p_w[j] for j in 1:4) for i in 1:4]
        [R_b_w[:, i] ~ quaternion_to_rotation_matrix(quaternion_multiply(Q_p_w, Q_b_p))[:, i] for i in 1:3] # https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation#Performance_comparisons

        D(ω_p[1]) ~ α_p[1]
        D(ω_p[2]) ~ α_p[2]
        D(ω_p[3]) ~ α_p[3]
        α_p[1] ~ (torque_p[1] + (s.I_kite[2] - s.I_kite[3]) * ω_p[2] * ω_p[3]) / s.I_kite[1] - 0.0s.orient_damping*ω_p[1]
        α_p[2] ~ (torque_p[2] + (s.I_kite[3] - s.I_kite[1]) * ω_p[3] * ω_p[1]) / s.I_kite[2] - 0.0s.orient_damping*ω_p[2]
        α_p[3] ~ (torque_p[3] + (s.I_kite[1] - s.I_kite[2]) * ω_p[1] * ω_p[2]) / s.I_kite[3] - 0.0s.orient_damping*ω_p[3]

        [D(kite_pos[i]) ~ kite_vel[i] for i in 1:3]
        [D(kite_vel[i]) ~ kite_acc[i] for i in 1:3]

        [pos[:, i]              .~ 0.0 for i in 1:3]
        [D.(pos[:, i])          .~ vel[:, i] for i in 4:s.i_A-1]
        D(trailing_edge_angle)   ~ trailing_edge_ω
        [vel[:, i]              .~ 0.0 for i in 1:3]
        [D.(vel[:, i])          .~ acc[:, i] for i in 4:s.i_A-1]
        D(trailing_edge_ω)       ~ trailing_edge_α
        D.(tether_length)       .~ tether_vel
        D.(tether_vel)          .~ tether_acc
    ]

    # Compute the masses and forces
    force_eqs = SizedArray{Tuple{3, s.i_C}, Symbolics.Equation}(undef)
    force_eqs[:, :] .= (force[:, :] .~ 0)
    
    seqs = []
    seqs            = scalar_eqs!(s, seqs, pos, vel, R_b_w, ω_p, ω_b, kite_pos, kite_vel, trailing_edge_angle, trailing_edge_ω, segment_length, mass_tether_particle, damping, c_spring, 
                        e_y, e_z, e_x, e_r_D, e_r_E, e_te_A, e_te_B, rho_kite, tether_length, tether_vel,
                        mass_per_meter, force, set_values)
    seqs, force_eqs = calc_kite_forces!(s, seqs, force_eqs, force, torque_p, R_b_w, kite_vel, kite_acc, ω_b, t, e_x, e_z, rho_kite, v_wind, trailing_edge_angle)
    seqs, force_eqs = calc_tether_forces!(s, seqs, force_eqs, t, force, pos, vel, segment_length, c_spring, damping, v_wind_gnd, norm1)
    
    if s.torque_control
        seqs = vcat(seqs, tether_acc .~ [calc_acc_torque(s.motors[i], tether_vel[i], norm(force[:, (i-1)%3+1]),
            set_values[i]) for i in 1:3])
    else
        seqs = vcat(seqs, tether_acc .~ [calc_acc_speed(s.motors[i], tether_vel[i], norm(force[:, (i-1)%3+1]),
            set_values[i]) for i in 1:3])
    end
    for i in 1:3
        seqs = vcat(seqs, vcat(force_eqs[:, i]))
        seqs = vcat(seqs, acc[:, i] .~ 0)
    end
    for i in 4:s.i_A-1
        seqs = vcat(seqs, vcat(force_eqs[:, i]))
        seqs = vcat(seqs, acc[:, i] .~ [0.0, 0.0, -G_EARTH] + (force[:, i] / mass_tether_particle[(i-1)%3+1]))
    end

    # torque = I * trailing_edge_α
    # trailing_edge_α = torque / (1/3 * (kite_mass/8) * kite_length_D^2)
    # torque = force[:, i] * kite_length_D
    # trailing_edge_α = force[:, i] * kite_length_D / (1/3 * (kite_mass/8) * kite_length_D^2)

    # 1. add all flap + spring + drag forces to flap_C point
    # 2. remove forces not in e_flap_c direction
    # 3. calculate acceleration from force flap c in e_flap_c direction

    te_I = (1/3 * (s.set.mass/8) * te_length^2)
    # -damping / I * ω = α_damping
    # solve for c: (c * (k*m/s^2) / (k*m^2)) * (m/s)=m/s^2 in wolframalpha
    # damping should be N*m*s
    rot_damping = 0.1s.damping * te_length

    seqs = [
        seqs
        vcat(force_eqs[:, s.i_A])
        vcat(force_eqs[:, s.i_B])
        vcat(force_eqs[:, s.i_C])
         # TODO: add turning drag instead of damping
        trailing_edge_α[1] ~ (force[:, s.i_A]) ⋅ e_te_A * te_length / te_I - (rot_damping[1] / te_I) * trailing_edge_ω[1]
        trailing_edge_α[2] ~ (force[:, s.i_B]) ⋅ e_te_B * te_length / te_I - (rot_damping[2] / te_I) * trailing_edge_ω[2]
    ]
    
    eqs = vcat(deqs, seqs)
    eqs = Symbolics.scalarize.(reduce(vcat, Symbolics.scalarize.(eqs)))

    discrete_events = [true => [Q_p_w[i] ~ normalize(Q_p_w)[i] for i in 1:4]]
    
    @named sys = ODESystem(eqs, t; discrete_events)
    return sys, collect(set_values)
end

"""
The distance/vel/acc of the kite cannot be measured directly, but the average acc of the kite distance is equal to the tether acc.
    So while the set_values input and wind is constant, distance_acc = tether_acc[3]. To accurately describe distance_acc when set_values
    or wind suddenly change, the distance_acc is calculated by combining the simulated distance_acc from timestep-1 to timestep and the
    current measured acc.
    distance_acc = 0.99 * sim_distance_acc + 0.01 * measured_distance_acc

Assume distance_acc = tether_acc[3] for convenience
"""
function model!(s::KPSQ; real=true)
    s.ϵ = 1e-6
    # pos, vel = init_pos_vel(s)
    init_pos!(s; α = 10.0)
    sys, inputs = create_sys!(s)
    # structural_simplify(sys, (inputs, []))
    (sys, _) = structural_simplify(sys, (inputs, []); fully_determined=false)
    s.simple_sys = sys

    u0map = [
        # tether pos --> expected_pos
        # tether vel --> s.tether_vel
        # tether acc --> 0

        # kite pos --> s.kite_pos
        # kite vel --> whatever
        # kite acc --> 0

        [sys.Q_p_w[i] => s.Q_p_w[i] for i in 1:4]
        # [sys.Q_vel[i] => 0 for i in 1:4]
        [sys.ω_p[i] => 0 for i in 1:3]
        # [sys.α_p[i] => 0 for i in 1:3]

        [sys.kite_pos[j] => s.kite_pos[j] for j in 1:3]
        [sys.kite_vel[i] => 0.0 for i in 1:3]
        # [sys.kite_acc[i] => 0.0 for i in 1:3]

        [sys.pos[j, i] => s.pos[j, i] for j in 1:3 for i in 4:s.i_A-1]
        [sys.vel[j, i] => 0 for j in 1:3 for i in 4:s.i_A-1]
        # [sys.acc[j, i] => 0 for j in 1:3 for i in 4:s.i_A-1]

        [sys.trailing_edge_angle[i] => 0.3 for i in 1:2]
        [sys.trailing_edge_ω[i] => 0 for i in 1:2]
        # [sys.trailing_edge_α[i] => 0 for i in 1:2]

        sys.gust_factor => 1.0

        [sys.tether_length[i] => s.measure.tether_length[i] for i in 1:3]
        [sys.tether_vel[j] => 0 for j in 1:3]
        # [sys.tether_acc[i] => 0 for i in 1:3]
    ]
    guesses = [
        # [sys.q[i] => s.q[i] for i in 1:4]
        # [sys.Q_vel[i] => 0 for i in 1:4]
        # [sys.ω_p[i] => 0 for i in 1:3]
        # [sys.α_p[i] => 0 for i in 1:3]

        # [sys.kite_pos[i] => s.kite_pos[i] for i in 1:3]
        # [sys.kite_vel[i] => 0 for i in 1:3]
        # [sys.kite_acc[i] => 0 for i in 1:3]

        # [sys.pos[j, i] => s.pos[j, i] for j in 1:3 for i in 4:s.i_C]
        # [sys.vel[j, i] => 0 for j in 1:3 for i in 4:s.i_C]
        # [sys.acc[j, i] => 0 for j in 1:3 for i in 4:s.i_C]

        # [sys.expected_pos[j, i] => s.pos[j, i] for j in 1:3 for i in 4:s.i_C]
        # [sys.expected_vel[j, i] => 0 for j in 1:3 for i in 4:s.i_C]

        # [sys.trailing_edge_angle[j] => 0 for j in 1:2]
        # [sys.trailing_edge_ω[j] => 0 for j in 1:2]
        # [sys.trailing_edge_α[j] => 0 for j in 1:2]

        # [sys.tether_length[j] => s.tether_lengths[j] for j in 1:3]
        # [sys.tether_vel[j] => 0 for j in 1:3]
        # [sys.tether_acc[j] => 0 for j in 1:3]

        # sys.gust_factor => 1.0natural
        ]
    p0map = [sys.set_values[j] => s.measure.winch_torque[j] for j in 1:3]
    # println("making prob")
    # @time prob = ModelingToolkit.InitializationProblem(sys, 0.0, u0map, p0map; guesses, fully_determined=true)

    # @show typeof(prob)
    # tol = 1e-3
    # # @time remake(prob; u0=u0map)
    # println("solving")
    # @time sol = NonlinearSolve.solve(prob)
    # println(sol.resid)
    # @show norm(sol.resid)
    # @show norm(sol[sys.q])
    
    # println("remaking init prob")
    # @time remake(prob; u0=[sys.gust_factor=>1.1])
    # @time sol = solve(prob; maxiters=10_000, abstol=tol, reltol=tol)
    # @time remake(prob; u0=[sys.gust_factor=>0.9])
    # @time sol = solve(prob; maxiters=10_000, abstol=tol, reltol=tol)
    
    dt = 1/s.set.sample_freq
    tspan   = (0.0, dt)
    # @time s.prob = ODEProblem(sys, u0map, tspan; guesses, fully_determined=real, check_length=real)
    # u0 = collect(unknowns(sys) .=> sol[unknowns(sys)])
    # p0 = collect(sys.set_values .=> sol[sys.set_values])
    s.prob = ODEProblem(sys, u0map, tspan, p0map)

    solver = QNDF(autodiff=false) # https://docs.sciml.ai/scimlbenchmarksoutput/stable/#results
    s.integrator = OrdinaryDiffEqCore.init(s.prob, solver; dt, abstol=s.set.abs_tol, reltol=s.set.rel_tol)
    return s.prob, sys, u0map
end