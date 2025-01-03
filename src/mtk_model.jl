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

function update_state!(s::KPS4_3L)
    pos = s.get_pos(s.integrator)
    [s.pos[i]          .= pos[:, i] for i in 1:s.i_A]
    s.tether_lengths   .= s.get_tether_lengths(s.integrator)
    s.vel_kite         .= s.get_vel_kite(s.integrator)
    calc_kite_ref_frame!(s, s.pos[s.i_E], s.pos[s.i_C], s.pos[s.i_D])
    # @assert all(abs.(s.trailing_edge_angle) .<= deg2rad(90))
    nothing
end

function convert_pos_vel(s::KPS4_3L, pos_, vel_)
    pos = Array{Union{Nothing, Float64}}(nothing, 3, s.i_A)
    vel = Array{Union{Nothing, Float64}}(nothing, 3, s.i_A)
    [pos[:,i] .= pos_[i] for i in 1:s.i_A-1]
    [vel[:,i] .= vel_[i] for i in 1:s.i_A-1]
    [pos[:,i] .= pos_[i] for i in s.i_B+1:s.i_A]
    [vel[:,i] .= vel_[i] for i in s.i_B+1:s.i_A]
    return pos, vel
end

function rotate_by_quaternion(vec, q)
    w, x, y, z = q[1], q[2], q[3], q[4]
    [1-2*(y^2+z^2)    2*(x*y-w*z)     2*(x*z+w*y);
     2*(x*y+w*z)      1-2*(x^2+z^2)   2*(y*z-w*x);
     2*(x*z-w*y)      2*(y*z+w*x)     1-2*(x^2+y^2)] * vec
end

function rotation_matrix_to_quaternion(R::Matrix)
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

function enu_to_ned(vec)
    return [-vec[1], vec[2], -vec[3]]
end

"""
Calculate the forces acting on the kite inertia particle.
"""
function calc_kite_forces!(s::KPS4_3L, seqs, force_eqs, force, torque_p, q, kite_pos, kite_vel, kite_acc, ω_p, t, e_x, e_z, rho, v_wind, trailing_edge_angle)
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
        total_torque_p(t)[1:3]
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
        circle_origin = (kite_pos + s.circle_center_t * e_z)
        seqs = [
            seqs
            aero_pos[:, i]  ~ kite_pos + rotate_by_quaternion(s.seg_cop_pos_b[:, i], q)
            e_r[:, i]       ~ normalize(circle_origin - aero_pos[:, i])
            seg_vel[:, i]   ~ kite_vel + rotate_by_quaternion(cross(ω_p, s.seg_cop_pos_p[:, i]), q)
            v_a[:, i]       ~ v_wind .- seg_vel[:, i]
            e_drift[:, i]   ~ (e_r[:, i] × -e_x)
            v_a_xr[:, i]    ~ v_a[:, i] .- (v_a[:, i] ⋅ e_drift[:, i]) .* e_drift[:, i]

            α < α_m_l ?
                seg_trailing_edge_angle[i] ~ trailing_edge_angle[1]  :
            α > α_m_r ?
                seg_trailing_edge_angle[i] ~ trailing_edge_angle[2] :
                seg_trailing_edge_angle[i] ~ ((trailing_edge_angle[2] - trailing_edge_angle[1]) / (α_m_r - α_m_l) * (α - α_m_l) + trailing_edge_angle[1])

            aoa[i]      ~ -asin2((v_a_xr[:, i] / norm(v_a_xr[:, i])) ⋅ e_r[:, i]) + deg2rad(s.set.alpha_zero)
            seg_cl[i]   ~ sym_interp(s.cl_interp, aoa[i], seg_trailing_edge_angle[i])
            seg_cd[i]   ~ sym_interp(s.cd_interp, aoa[i], seg_trailing_edge_angle[i])

            seg_L[:, i] ~ 0.5 * rho * (norm(v_a_xr[:, i]))^2 * s.set.radius * dα * kite_length * seg_cl[i] * 
                                ((v_a_xr[:, i] × e_drift[:, i]) / norm(v_a_xr[:, i] × e_drift[:, i]))
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
            seg_aero_torque_p[:, i] ~ rotate_by_quaternion(s.seg_cop_pos_p[:, i], q) × (s.R_b_p * seg_aero_force[:, i])
            seg_gravity_torque_p[:, i] ~ rotate_by_quaternion(s.seg_com_pos_p[:, i], q) × (s.R_b_p * seg_g_force[:, i])
        ]

        # if i <= n
        #     [f_te_c_eq[j] = (F_te_C[j] ~ f_te_c_eq[j].rhs + seg_te_force[j, i]) for j in 1:3]
        # else 
        #     [f_te_d_eq[j] = (F_te_D[j] ~ f_te_d_eq[j].rhs + seg_te_force[j, i]) for j in 1:3]
        # end
    end

    seqs = [
        seqs
        total_kite_force ~ [sum(seg_aero_force[i, :]) + sum(seg_g_force[i, :]) + force[i, s.i_C] for i in 1:3]
        kite_acc ~ total_kite_force / s.set.mass
        C_torque_p ~ rotate_by_quaternion(s.pos_C_p, q) × (s.R_b_p * force[:, s.i_C])
        total_torque_p ~ [sum(seg_aero_torque_p[i, :]) + sum(seg_gravity_torque_p[i, :]) + C_torque_p[i] for i in 1:3]
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
    calc_particle_forces!(s::KPS4_3L, seqs, force_eqs, force, pos1, pos2, vel1, vel2, length, c_spring, damping, 
                          rho, i, l_0, k, c, segment, rel_vel, av_vel, norm1, unit_vector, k1, k2, c1, spring_vel,
                          spring_force, v_apparent, v_wind_tether, area, v_app_perp, half_drag_force)

Calculate the drag force and spring force of the tether segment, defined by the parameters pos1, pos2, vel1 and vel2
and distribute it equally on the two particles, that are attached to the segment.
The result is stored in the array s.forces. 
"""
function calc_particle_forces!(s::KPS4_3L, seqs, force_eqs, force, p1, p2, pos1, pos2, vel1, vel2, length, c_spring, 
    damping, rho, i, l_0, k, c, segment, rel_vel, av_vel, norm1, unit_vector, k1, k2, c1, c2, spring_vel, perp_vel,
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
        k1           ~ 1.0 * k # compression stiffness kite segments
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
                ((k  * (l_0 - norm1) - c * spring_vel) * unit_vector[j] - c2 * perp_vel[j]) * (1 + smooth_sign_ϵ(norm1 - l_0; s.ϵ)) / 2 +
                ((-c * spring_vel) * unit_vector[j] - c2 * perp_vel[j]) * (1 - smooth_sign_ϵ(norm1 - l_0; s.ϵ)) / 2
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
@inline function calc_tether_forces!(s::KPS4_3L, seqs, force_eqs, t, force, pos, vel, length, c_spring, damping, v_wind_gnd, norm1)
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
        k1(t)[1:s.i_C-3]
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
                          rel_vel[:, i], av_vel[:, i], norm1[i], unit_vector[:, i], k1[i], k2[i], c1[i], c2[i], spring_vel[i],
                          perp_vel[:, i], spring_force[:, i], v_apparent[:, i], v_wind_tether[:, i], area[i], v_app_perp[:, i], 
                          half_drag_force[:, i])
    end

    return seqs, force_eqs
end

function expected_pos_vel(s, seqs, pos, kite_pos, kite_vel, tether_vel, tether_length, tether_force, norm1)
    @variables begin
        stretched_tether_length(t)[1:3]
        expected_pos(t)[1:3, 1:s.i_C]
        tether_move_vel(t)[1:3, 1:s.i_C]
        tether_kite_vel(t)[1:3, 1:s.i_C]
        expected_vel(t)[1:3, 1:s.i_C]
    end
    seqs = [
        seqs
        [stretched_tether_length[i] ~ sum([norm1[j] for j in range(i, step=3, length=s.set.segments)]) for i in 1:3]
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

function scalar_eqs!(s, seqs, pos, vel, q, ω_p, ω_b, kite_pos, kite_vel, trailing_edge_angle, trailing_edge_ω, segment_length, mass_tether_particle, damping, c_spring, 
        e_y, e_z, e_x, e_r_D, e_r_E, e_te_C, e_te_D, rho_kite, damping_coeff, tether_length, tether_vel,
        mass_per_meter, force, set_values, norm1)
    
    te_length = s.kite_length_D/4
    # last tether point vel without rotational trailing edge vel
    vel_kite_A = kite_vel + rotate_by_quaternion(cross(ω_p, s.pos_A_p), q)
    vel_kite_B = kite_vel + rotate_by_quaternion(cross(ω_p, s.pos_B_p), q)
    vel_kite_C = kite_vel + rotate_by_quaternion(cross(ω_p, s.pos_C_p), q)
    seqs = [
        seqs
        ω_b ~ s.R_b_p' * ω_p
        
        # last tether point pos
        pos[:, s.i_C]    ~ kite_pos - e_z * s.C_t
        pos[:, s.i_A]    ~ pos[:, s.i_C] + e_y * s.pos_D_b[2] + e_x * te_length * cos(trailing_edge_angle[1]) + 
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
        e_y     ~ rotate_by_quaternion([0, 1, 0], q)
        e_z     ~ rotate_by_quaternion([0, 0, 1], q)
        e_x     ~ rotate_by_quaternion([1, 0, 0], q)
        e_r_D   ~ normalize(rotate_by_quaternion(-s.pos_D_b, q))
        e_r_E   ~ normalize(rotate_by_quaternion(-s.pos_E_b, q))
        e_te_C  ~ -e_x * sin(trailing_edge_angle[1]) + e_r_D * cos(trailing_edge_angle[1])
        e_te_D  ~ -e_x * sin(trailing_edge_angle[2]) + e_r_E * cos(trailing_edge_angle[2])
        rho_kite        ~ calc_rho(s.am, pos[3,s.i_A])
        damping_coeff   ~ max(1.0 - t, 0.0) * s.damping_coeff
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
    # println(size([(vel[:, i] ⋅ normalize(pos[:, i]) for i in s.i_A: s.i_E)]))
    # println([kite_vel[i] ~ vel[:, j] ⋅ normalize(pos[:, j]) for (i, j) in enumerate(s.i_A: s.i_E)])
    seqs = expected_pos_vel(s, seqs, pos, kite_pos, kite_vel, tether_vel, tether_length, tether_force, norm1)
    seqs = [
        seqs
        tether_force     ~ [norm(force[:, i]) for i in 1:3]
        # orientation     ~ orient_euler(s; one_point=false)
        # upwind_dir      ~ upwind_dir(v_wind_gnd)
        # heading         ~ calc_heading(orientation, elevation, azimuth; upwind_dir)
        heading_y       ~ calc_heading_y(-e_x)
        power_angle         ~ (trailing_edge_angle[1] + trailing_edge_angle[2]) / 2
        power_vel           ~ (trailing_edge_ω[1] + trailing_edge_ω[2]) / 2
        steering_angle      ~ trailing_edge_angle[2] - trailing_edge_angle[1]
        steering_vel        ~ trailing_edge_ω[2] - trailing_edge_ω[1]
        tether_diff         ~ tether_length[2] - tether_length[1]
        tether_diff_vel     ~ tether_vel[2] - tether_vel[1]
        set_diff            ~ set_values[2] - set_values[1]
        # D(set_values)       ~ [0, 0, 0]

        distance_vel        ~ kite_vel ⋅ normalize(kite_pos - pos[:, 3])
        distance_acc        ~ kite_acc ⋅ normalize(kite_pos - pos[:, 3])
        elevation           ~ atan(kite_pos[3] / kite_pos[1])
        elevation_vel       ~ D(elevation)
        azimuth             ~ -atan(kite_pos[2] / kite_pos[1])
        azimuth_vel         ~ D(azimuth)
        x_acc               ~ kite_acc ⋅ e_x
        y_acc               ~ kite_acc ⋅ e_y
        left_diff           ~ tether_length[1] - tether_length[3]
        right_diff          ~ tether_length[2] - tether_length[3]
        [distance[i]        ~ norm(pos[:, i]) for i in 1:s.i_C]
    ]
    return seqs
end

function create_sys!(s::KPS4_3L)
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
        q(t)[1:4] # quaternion orientation of the kite
        ω_p(t)[1:3] # turn rate in principal frame
        ω_b(t)[1:3] # turn rate in body frame
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
        damping_coeff(t)
        c_spring(t)[1:3]
        e_x(t)[1:3]
        e_y(t)[1:3]
        e_z(t)[1:3]
        e_r_D(t)[1:3]
        e_r_E(t)[1:3]
        e_te_C(t)[1:3]
        e_te_D(t)[1:3]
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

    deqs = [
        deqs
        [D(q[i]) ~ 0.5 * sum(Ω[i, j] * q[j] for j in 1:4) for i in 1:4]
        D(ω_p[1]) ~ (torque_p[1] + (s.I_kite[2] - s.I_kite[3]) * ω_p[2] * ω_p[3]) / s.I_kite[1] - s.orient_damping*ω_p[1]
        D(ω_p[2]) ~ (torque_p[2] + (s.I_kite[3] - s.I_kite[1]) * ω_p[3] * ω_p[1]) / s.I_kite[2] - s.orient_damping*ω_p[2]
        D(ω_p[3]) ~ (torque_p[3] + (s.I_kite[1] - s.I_kite[2]) * ω_p[1] * ω_p[2]) / s.I_kite[3] - s.orient_damping*ω_p[3]
        [D(kite_pos[i]) ~ kite_vel[i] for i in 1:3]
        [D(kite_vel[i]) ~ kite_acc[i] for i in 1:3]

        [pos[:, i] .~ 0.0 for i in 1:3]
        [D.(pos[:, i]) .~ vel[:, i] for i in 4:s.i_A-1]
        D(trailing_edge_angle)   ~ trailing_edge_ω
        [vel[:, i] .~ 0.0 for i in 1:3]
        [D.(vel[:, i]) .~ acc[:, i] for i in 4:s.i_A-1]
        D(trailing_edge_ω)   ~ trailing_edge_α
        D.(tether_length) .~ tether_vel
        D.(tether_vel) .~ tether_acc
    ]

    # Compute the masses and forces
    force_eqs = SizedArray{Tuple{3, s.i_C}, Symbolics.Equation}(undef)
    force_eqs[:, :] .= (force[:, :] .~ 0)
    
    seqs = []
    seqs            = scalar_eqs!(s, seqs, pos, vel, q, ω_p, ω_b, kite_pos, kite_vel, trailing_edge_angle, trailing_edge_ω, segment_length, mass_tether_particle, damping, c_spring, 
                        e_y, e_z, e_x, e_r_D, e_r_E, e_te_C, e_te_D, rho_kite, damping_coeff, tether_length, tether_vel,
                        mass_per_meter, force, set_values, norm1)
    seqs, force_eqs = calc_kite_forces!(s, seqs, force_eqs, force, torque_p, q, kite_pos, kite_vel, kite_acc, ω_p, t, e_x, e_z, rho_kite, v_wind, trailing_edge_angle)
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
        seqs = vcat(seqs, acc[:, i] .~ [0.0; 0.0; -G_EARTH] .+ (force[:, i] ./ mass_tether_particle[(i-1)%3+1]) .- damping_coeff * vel[:, i])
    end

    # torque = I * trailing_edge_α
    # trailing_edge_α = torque / (1/3 * (kite_mass/8) * kite_length_D^2)
    # torque = force[:, i] * kite_length_D
    # trailing_edge_α = force[:, i] * kite_length_D / (1/3 * (kite_mass/8) * kite_length_D^2)

    # 1. add all flap + spring + drag forces to flap_C point
    # 2. remove forces not in e_flap_c direction
    # 3. calculate acceleration from force flap c in e_flap_c direction

    te_I = (1/3 * (s.set.mass/8) * te_length^2)
    critical_damping = 2.0 * sqrt(te_I)
    @show critical_damping
    seqs = [
        seqs
        vcat(force_eqs[:, s.i_A])
        vcat(force_eqs[:, s.i_B])
         # TODO: add turning drag instead of damping
        trailing_edge_α[1] ~ (force[:, s.i_A]) ⋅ e_te_C * te_length / te_I - critical_damping * trailing_edge_ω[1]
            # (damping_coeff*200) * trailing_edge_ω[1]
        trailing_edge_α[2] ~ (force[:, s.i_B]) ⋅ e_te_D * te_length / te_I - critical_damping * trailing_edge_ω[2]
            # (damping_coeff*200) * trailing_edge_ω[1]
    ]

    eqs = vcat(deqs, seqs)
    eqs = Symbolics.scalarize.(reduce(vcat, Symbolics.scalarize.(eqs)))

    @named sys = ODESystem(eqs, t)
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
function model!(s::KPS4_3L; real=true)
    # pos, vel = init_pos_vel(s)
    init_pos!(s; α = 10.0)
    sys, inputs = create_sys!(s)
    (sys, _) = structural_simplify(sys, (inputs, []); fully_determined=false)
    s.simple_sys = sys

    u0map = []
    @show s.pos
    if real
        u0map = [
            # tether pos --> expected_pos
            # tether vel --> s.tether_vel
            # tether acc --> 0

            # kite pos --> s.kite_pos
            # kite vel --> whatever
            # kite acc --> 0

            [sys.kite_pos[j] => s.kite_pos[j] for j in 1:3]
            [sys.kite_vel[j] => 0.0 for j in 1:3]

            # [sys.pos[j, s.i_A] => s.pos[s.i_A][j] for j in 1:3]
            # [sys.vel[j, s.i_A] => [-1, 0, 0][j] for j in 1:3]

            [sys.vel[j, i] => 0 for j in 1:3 for i in 4:s.i_A-1]
            [sys.acc[j, i] => 0 for j in 1:3 for i in 4:s.i_A-1]

            [sys.trailing_edge_ω[j] => 0 for j in 1:2]
            [sys.trailing_edge_α[j] => 0 for j in 1:2]

            sys.gust_factor => 1.0

            [sys.tether_length[j] => s.measure.tether_length[j] for j in 1:3]
            # [sys.tether_vel[j] => [1.46, 1.46, 1.46][j] for j in 1:3]
            [sys.tether_acc[j] => 0 for j in 1:3]
        ]
    else
        u0map = [
            # sys.distance[s.i_A] => norm(s.pos[s.i_A])
            [sys.pos[j, s.i_A] => s.pos[s.i_A][j] for j in 1:3]
            [sys.vel[j, i] => norm(s.pos[i]) / norm(s.pos[s.i_A]) * sys.vel[j, s.i_A] 
                for j in 1:3 for i in vcat(4:s.i_A-1, s.i_B+1:s.i_A-1)]
            [sys.acc[j, i] => 0 for j in 1:3 for i in vcat(4:s.i_A-1, s.i_B+1:s.i_A)]

            [sys.trailing_edge_ω[j] => 0 for j in 1:2]
            [sys.trailing_edge_α[j] => 0 for j in 1:2]

            sys.gust_factor => 1.0

            [sys.tether_length[j] => s.measure.tether_length[j] for j in 1:3]
            # [sys.tether_vel[j] => [1.46, 1.46, 1.46][j] for j in 1:3]
            [sys.tether_acc[j] => 0 for j in 1:3]
        ]
    end
    guesses = [
        [sys.kite_pos[j] => s.kite_pos[j] for j in 1:3]
        [sys.kite_vel[j] => 0 for j in 1:3]
        [sys.kite_acc[j] => 0 for j in 1:3]
        [sys.pos[j, i] => s.pos[j, i] for j in 1:3 for i in 4:s.i_C]
        [sys.expected_pos[j, i] => s.pos[j, i] for j in 1:3 for i in 4:s.i_C]
        [sys.vel[j, i] => 0 for j in 1:3 for i in 4:s.i_C]
        [sys.expected_vel[j, i] => 0 for j in 1:3 for i in 4:s.i_C]
        [sys.acc[j, i] => 0 for j in 1:3 for i in 4:s.i_C]

        [sys.trailing_edge_angle[j] => 0 for j in 1:2]
        [sys.trailing_edge_ω[j] => 0 for j in 1:2]

        [sys.tether_length[j] => s.tether_lengths[j] for j in 1:3]
        [sys.tether_vel[j] => 0 for j in 1:3]

        sys.gust_factor => 1.0
        [sys.set_values[j] => s.measure.winch_torque[j] for j in 1:3]
        sys.set_diff => s.measure.winch_torque[2] - s.measure.winch_torque[1]
    ]
    @time prob = ModelingToolkit.InitializationProblem(sys, 0.0, u0map; guesses, fully_determined=true)
    tol = 1e-3
    # @time remake(prob; u0=u0map)
    # @time sol = solve(prob, RobustMultiNewton(); maxiters=10_000, abstol=tol, reltol=tol)
    @time sol = solve(prob; maxiters=20_000, abstol=tol, reltol=tol)
    @time sol = solve(prob; maxiters=20_000, abstol=tol, reltol=tol)
    @time sol = solve(prob; maxiters=20_000, abstol=tol, reltol=tol)
    # println("remaking init prob")
    # @time remake(prob; u0=[sys.gust_factor=>1.1])
    # @time sol = solve(prob; maxiters=10_000, abstol=tol, reltol=tol)
    # @time remake(prob; u0=[sys.gust_factor=>0.9])
    # @time sol = solve(prob; maxiters=10_000, abstol=tol, reltol=tol)

    dt = 1/s.set.sample_freq
    tspan   = (0.0, dt)
    # @time s.prob = ODEProblem(sys, u0map, tspan; guesses, fully_determined=real, check_length=real)
    u0 = collect(unknowns(sys) .=> sol[unknowns(sys)])
    p0 = collect(sys.set_values .=> sol[sys.set_values])
    s.prob = ODEProblem(sys, u0, tspan, p0)
    solver = QNDF(autodiff=false) # https://docs.sciml.ai/scimlbenchmarksoutput/stable/#results
    s.integrator = OrdinaryDiffEqCore.init(s.prob, solver; dt, abstol=s.set.abs_tol, reltol=s.set.rel_tol, save_on=false)
    return prob, sol, sys, u0map
end