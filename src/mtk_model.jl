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

function sym_interp(interp::Function, aoa, flap_angle)
    return interp(rad2deg(aoa), rad2deg(flap_angle-aoa))
end
@register_symbolic sym_interp(interp::Function, aoa, flap_angle)

function normalize(vec)
    return vec / norm(vec)
end
@register_symbolic normalize(vec)

function update_state!(s)
    pos = s.get_pos(s.integrator)
    [s.pos[i]          .= pos[:, i] for i in 1:s.num_A]
    s.tether_lengths   .= s.get_tether_lengths(s.integrator)
    s.vel_kite         .= s.get_vel_kite(s.integrator)
    calc_kite_ref_frame!(s, s.pos[s.num_E], s.pos[s.num_C], s.pos[s.num_D])
    # @assert all(abs.(s.flap_angle) .<= deg2rad(90))
    nothing
end

function convert_pos_vel(s::KPS4_3L, pos_, vel_)
    pos = Array{Union{Nothing, Float64}}(nothing, 3, s.num_A)
    vel = Array{Union{Nothing, Float64}}(nothing, 3, s.num_A)
    [pos[:,i] .= pos_[i] for i in 1:s.num_flap_C-1]
    [vel[:,i] .= vel_[i] for i in 1:s.num_flap_C-1]
    [pos[:,i] .= pos_[i] for i in s.num_flap_D+1:s.num_A]
    [vel[:,i] .= vel_[i] for i in s.num_flap_D+1:s.num_A]
    return pos, vel
end

"""
    calc_aero_forces!(s::KPS4_3L, seqs, force_eqs, force, pos, vel, t, e_x, e_y, e_z, rho)

Calculates the aerodynamic forces acting on the kite particles.

Parameters:
- pos:              vector of the particle positions
- vel:              vector of the particle velocities
- rho:              air density [kg/m^3]

Updates the vector s.forces of the first parameter.
"""
function calc_aero_forces!(s::KPS4_3L, seqs, force_eqs, force, pos, vel, t, e_x, e_y, e_z, E_c, rho, v_wind, flap_angle)
    n = s.set.aero_surfaces
    @variables begin
        v_cx(t)[1:3]
        v_dx(t)[1:3]
        v_dy(t)[1:3]
        v_dz(t)[1:3]
        v_cy(t)[1:3]
        v_cz(t)[1:3]
        y_lc(t)
        y_ld(t)
    end

    seqs = [
        seqs
        v_cx    ~ (vel[:, s.num_C] ⋅ e_x) * e_x
        v_dx    ~ (vel[:, s.num_D] ⋅ e_x) * e_x
        v_dy    ~ (vel[:, s.num_D] ⋅ e_y) * e_y
        v_dz    ~ (vel[:, s.num_D] ⋅ e_z) * e_z
        v_cy    ~ (vel[:, s.num_C] ⋅ e_y) * e_y
        v_cz    ~ (vel[:, s.num_C] ⋅ e_z) * e_z
        y_lc    ~  norm(pos[:, s.num_C] - 0.5 * (pos[:, s.num_C] + pos[:, s.num_D]))
        y_ld    ~ -norm(pos[:, s.num_D] - 0.5 * (pos[:, s.num_C] + pos[:, s.num_D]))
    ]

    # integrating loop variables, iterating over 2n segments
    @variables begin
        F(t)[1:3, 1:2n]
        e_r(t)[1:3, 1:2n]
        e_te(t)[1:3, 1:2n] # clockwise trailing edge of flap vector
        y_l(t)[1:2n]
        v_kite(t)[1:3, 1:2n]
        v_a(t)[1:3, 1:2n]
        e_drift(t)[1:3, 1:2n]
        v_a_xr(t)[1:3, 1:2n]
        aoa(t)[1:n*2]
        cl_seg(t)[1:n*2]
        cd_seg(t)[1:n*2]
        L_seg(t)[1:3, 1:2n]
        D_seg(t)[1:3, 1:2n]
        F_te_seg(t)[1:3, 1:2n]
        seg_flap_angle(t)[1:2n] # flap angle relative to -e_x, e_z
        ram_force(t)[1:2n]
        te_force(t)[1:2n]
        L_C(t)[1:3]
        L_D(t)[1:3]
        D_C(t)[1:3]
        D_D(t)[1:3]
        F_te_C(t)[1:3] # trailing edge clockwise force
        F_te_D(t)[1:3]
    end
    l_c_eq = SizedArray{Tuple{3}, Symbolics.Equation}(collect(L_C .~ 0))
    l_d_eq = SizedArray{Tuple{3}, Symbolics.Equation}(collect(L_D .~ 0))
    d_c_eq = SizedArray{Tuple{3}, Symbolics.Equation}(collect(D_C .~ 0))
    d_d_eq = SizedArray{Tuple{3}, Symbolics.Equation}(collect(D_D .~ 0))
    f_te_c_eq = SizedArray{Tuple{3}, Symbolics.Equation}(collect(F_te_C .~ 0))
    f_te_d_eq = SizedArray{Tuple{3}, Symbolics.Equation}(collect(F_te_D .~ 0))
    kite_length = zero(SimFloat)
    α           = zero(SimFloat)
    α_0         = zero(SimFloat)
    α_middle    = zero(SimFloat)
    dα          = zero(SimFloat)
    α_0         = π/2 - s.set.width/2/s.set.radius
    α_middle    = π/2
    dα          = (α_middle - α_0) / n
    for i in 1:n*2
        if i <= n
            α = α_0 + -dα/2 + i * dα
            kite_length = s.set.tip_length + (s.set.middle_length-s.set.tip_length) * (α - α_0) / (π/2 - α_0) # TODO: kite length gets less with flap turning
        else
            α = pi - (α_0 + -dα/2 + (i-n) * dα)
            kite_length = s.set.tip_length + (s.set.middle_length-s.set.tip_length) * (π - α_0 - α) / (π/2 - α_0)
        end
        seg_flap_height = kite_length * s.set.flap_height
        seqs = [
            seqs
            F[:, i]          ~ E_c + e_y * cos(α) * s.set.radius - e_z * sin(α) * s.set.radius
            e_r[:, i]        ~ (E_c - F[:, i]) / norm(E_c - F[:, i])
            y_l[i]           ~ cos(α) * s.set.radius
            α < π/2 ?
                v_kite[:, i] ~ ((v_cx - v_dx) / (y_lc - y_ld) * (y_l[i] - y_ld) + v_dx) + v_cy + v_cz :
                v_kite[:, i] ~ ((v_cx - v_dx) / (y_lc - y_ld) * (y_l[i] - y_ld) + v_dx) + v_dy + v_dz
            v_a[:, i]         ~ v_wind .- v_kite[:, i]
            e_drift[:, i]    ~ (e_r[:, i] × e_x)
            v_a_xr[:, i]     ~ v_a[:, i] .- (v_a[:, i] ⋅ e_drift[:, i]) .* e_drift[:, i]

            α < s.α_l ?
                seg_flap_angle[i]    ~ flap_angle[1]  :
            α > s.α_r ?
                seg_flap_angle[i]    ~ flap_angle[2] :
                seg_flap_angle[i]    ~ ((flap_angle[2] - flap_angle[1]) / (s.α_r - s.α_l) * (α - s.α_l) + (flap_angle[1]))

            aoa[i]      ~ -asin2((v_a_xr[:, i] / norm(v_a_xr[:, i])) ⋅ e_r[:, i]) + deg2rad(s.set.alpha_zero)
            cl_seg[i]   ~ sym_interp(s.cl_interp, aoa[i], seg_flap_angle[i])
            cd_seg[i]   ~ sym_interp(s.cd_interp, aoa[i], seg_flap_angle[i])

            L_seg[:, i] ~ 0.5 * rho * (norm(v_a_xr[:, i]))^2 * s.set.radius * dα * kite_length * cl_seg[i] * 
                                ((v_a_xr[:, i] × e_drift[:, i]) / norm(v_a_xr[:, i] × e_drift[:, i]))
            D_seg[:, i] ~ 0.5 * rho * norm(v_a_xr[:, i]) * s.set.radius * dα * kite_length * cd_seg[i] *
                                v_a_xr[:, i]


            e_te[:, i] ~ e_x * sin(seg_flap_angle[i]) + e_r[:, i] * cos(seg_flap_angle[i])
            ram_force[i] ~ smooth_sign_ϵ(deg2rad(s.set.alpha_zero) - seg_flap_angle[i]; s.ϵ) *
                        rho * norm(v_a[:, i])^2 * seg_flap_height * s.set.radius * dα * (seg_flap_height/2) / (kite_length/4)
            te_force[i] ~ 0.5 * rho * (norm(v_a_xr[:, i]))^2 * s.set.radius * dα * kite_length * 
                                sym_interp(s.c_te_interp, aoa[i], seg_flap_angle[i])
            F_te_seg[:, i] ~ (ram_force[i] + te_force[i]) * e_te[:, i]
        ]

        # TODO: correct for extra torque in wingtips (add to c substract from d)
        # TODO: use SymbolicNumericIntegration.jl
        if i <= n
            [l_c_eq[j] = (L_C[j] ~ l_c_eq[j].rhs + L_seg[j, i]) for j in 1:3]
            [d_c_eq[j] = (D_C[j] ~ d_c_eq[j].rhs + D_seg[j, i]) for j in 1:3]
            [f_te_c_eq[j] = (F_te_C[j] ~ f_te_c_eq[j].rhs + F_te_seg[j, i]) for j in 1:3]
        else 
            [l_d_eq[j] = (L_D[j] ~ l_d_eq[j].rhs + L_seg[j, i]) for j in 1:3]
            [d_d_eq[j] = (D_D[j] ~ d_d_eq[j].rhs + D_seg[j, i]) for j in 1:3]
            [f_te_d_eq[j] = (F_te_D[j] ~ f_te_d_eq[j].rhs + F_te_seg[j, i]) for j in 1:3]
        end
    end

    seqs = [
        seqs
        l_c_eq
        d_c_eq
        l_d_eq
        d_d_eq
        f_te_c_eq
        f_te_d_eq
    ]

    # longtitudinal force
    # F_inside_flap = P * A
    # F_inside_flap = rho * norm(v)^2 * flap_height * width
    # F_inside_flap = rho * norm(v)^2 * flap_height * radius * dα
    # F_trailing_edge = -F_inside_flap * (flap_height/2) / (kite_length/4) if flap_angle > 0 clockwise force
    # F_trailing_edge = F_inside_flap * (flap_height/2) / (kite_length/4) if flap_angle < 0 clockwise force
    # dF_te_dα = rho * norm(v)^2 * flap_height * radius
    # flap_height = height_middle * kite_length / middle_length
    
    # TODO: check if the right forces are added
    force_eqs[:,s.num_C]   .= (force[:,s.num_C]   .~ L_C + D_C) 
    force_eqs[:,s.num_D]   .= (force[:,s.num_D]   .~ L_D + D_D)
    force_eqs[:,s.num_flap_C] .= (force[:,s.num_flap_C] .~ F_te_C + [0.0, 0.0, -G_EARTH*(s.set.mass/8)])
    force_eqs[:,s.num_flap_D] .= (force[:,s.num_flap_D] .~ F_te_D + [0.0, 0.0, -G_EARTH*(s.set.mass/8)])
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
function calc_particle_forces!(s::KPS4_3L, seqs, force_eqs, force, pos1, pos2, vel1, vel2, length, c_spring, 
    damping, rho, i, l_0, k, c, segment, rel_vel, av_vel, norm1, unit_vector, k1, k2, c1, c2, spring_vel, perp_vel,
            spring_force, v_apparent, v_wind_tether, area, v_app_perp, half_drag_force)
    d_tether = s.set.d_tether/1000.0
    seqs = [
        seqs
        i <= s.set.segments*3 ? l_0 ~ length[(i-1) % 3 + 1] : l_0 ~ s.springs[i].length # Unstressed length
        i <= s.set.segments*3 ? k   ~ c_spring[(i-1) % 3 + 1] :
                                k   ~ s.springs[i].c_spring        # Spring constant
        i <= s.set.segments*3 ? c   ~ damping[(i-1) % 3 + 1] : c ~ s.springs[i].damping # Damping coefficient    
        segment     .~ pos1 - pos2 # TODO: all segments have same length and tension
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

    if i >= Base.length(s.springs) - KITE_SPRINGS_3L + 1  # kite springs
        for j in 1:3
            seqs = [
                seqs
                spring_force[j] ~ 
                    (k  * (l_0 - norm1) - c1 * spring_vel) * unit_vector[j] * (1 + smooth_sign_ϵ(norm1 - l_0; s.ϵ)) / 2 +
                    (k1 * (l_0 - norm1) -  c * spring_vel) * unit_vector[j] * (1 - smooth_sign_ϵ(norm1 - l_0; s.ϵ)) / 2
            ]
        end
    else
        for j in 1:3
            seqs = [
                seqs
                spring_force[j] ~
                    ((k  * (l_0 - norm1) - c * spring_vel) * unit_vector[j] - c2 * perp_vel[j]) * (1 + smooth_sign_ϵ(norm1 - l_0; s.ϵ)) / 2 +
                    ((-c * spring_vel) * unit_vector[j] - c2 * perp_vel[j]) * (1 - smooth_sign_ϵ(norm1 - l_0; s.ϵ)) / 2
            ]
        end
    end
    seqs = [
        seqs
        v_apparent       ~ v_wind_tether - av_vel
        i >= s.num_flap_C ?
            area             ~ norm1 * d_tether * 10 : # 10 is the number of parallel lines in the bridle system
            area             ~ norm1 * d_tether * (1 + (i%3 == 0)) # double area for middle tether
        v_app_perp       ~ v_apparent - (v_apparent ⋅ unit_vector) * unit_vector
        half_drag_force .~ (0.25 * rho * s.set.cd_tether * norm(v_app_perp) * area) .* v_app_perp
    ]

    for j in 1:3
        force_eqs[j, s.springs[i].p1] = 
            (force[j, s.springs[i].p1] ~ force_eqs[j, s.springs[i].p1].rhs + (half_drag_force[j] + spring_force[j]))
        force_eqs[j, s.springs[i].p2] = 
            (force[j, s.springs[i].p2] ~ force_eqs[j, s.springs[i].p2].rhs + (half_drag_force[j] - spring_force[j]))
    end
    
    return seqs, force_eqs
end



"""
    inner_loop!(s::KPS4_3L, seqs, force_eqs, t, force, pos, vel, length, c_spring, damping, v_wind_gnd)

Calculate the forces, acting on all particles.

Output:length
- s.forces
- s.v_wind_tether
"""
@inline function inner_loop!(s::KPS4_3L, seqs, force_eqs, t, force, pos, vel, length, c_spring, damping, v_wind_gnd)
    @variables begin
        height(t)[eachindex(s.springs)]
        rho(t)[eachindex(s.springs)]
        v_wind_tether(t)[1:3, eachindex(s.springs)]

        l_0(t)[eachindex(s.springs)]
        k(t)[eachindex(s.springs)]
        c(t)[eachindex(s.springs)]
        segment(t)[1:3, eachindex(s.springs)]
        rel_vel(t)[1:3, eachindex(s.springs)]
        av_vel(t)[1:3, eachindex(s.springs)] 
        norm1(t)[eachindex(s.springs)]
        unit_vector(t)[1:3, eachindex(s.springs)]
        k1(t)[eachindex(s.springs)]
        k2(t)[eachindex(s.springs)]
        c1(t)[eachindex(s.springs)]
        c2(t)[eachindex(s.springs)]
        spring_vel(t)[eachindex(s.springs)]
        perp_vel(t)[1:3, eachindex(s.springs)]
        spring_force(t)[1:3, eachindex(s.springs)]
        v_apparent(t)[1:3, eachindex(s.springs)]
        area(t)[eachindex(s.springs)]
        v_app_perp(t)[1:3, eachindex(s.springs)]
        half_drag_force(t)[1:3, eachindex(s.springs)]
        gust_factor(t)
    end
    
    seqs = [seqs; D(gust_factor) ~ 0]
    for i in eachindex(s.springs)
        p1 = s.springs[i].p1  # First point nr.
        p2 = s.springs[i].p2
        seqs = [
            seqs
            height[i]           ~ max(0.0, 0.5 * (pos[:, p1][3] + pos[:, p2][3]))
            rho[i]              ~ calc_rho(s.am, height[i])
            v_wind_tether[:, i] ~ AtmosphericModels.calc_wind_factor(s.am, height[i], s.set.profile_law) * v_wind_gnd * abs(gust_factor)
        ]

        seqs, force_eqs = calc_particle_forces!(s, seqs, force_eqs, force, pos[:, p1], pos[:, p2], vel[:, p1], 
                          vel[:, p2], length, c_spring, damping, rho[i], i, l_0[i], k[i], c[i], segment[:, i], 
                          rel_vel[:, i], av_vel[:, i], norm1[i], unit_vector[:, i], k1[i], k2[i], c1[i], c2[i], spring_vel[i],
                          perp_vel[:, i], spring_force[:, i], v_apparent[:, i], v_wind_tether[:, i], area[i], v_app_perp[:, i], 
                          half_drag_force[:, i])
    end

    return seqs, force_eqs
end

function scalar_eqs(s, seqs, pos, vel, acc, flap_angle, flap_vel, flap_acc, segment_length, mass_tether_particle, damping, c_spring, 
        P_c, e_y, e_z, e_x, e_r_C, e_r_D, e_te_C, e_te_D, E_c, rho_kite, damping_coeff, v_wind_gnd, tether_length, tether_vel,
        mass_per_meter, force, set_values)
    flap_length = s.kite_length_C/4
    seqs = [
        seqs
        pos[:, s.num_flap_C]    ~ pos[:, s.num_C] - e_x * flap_length * cos(flap_angle[1]) + e_r_C * flap_length * sin(flap_angle[1]) + 
            e_z * (sin(s.α_C) * s.set.radius + (s.set.bridle_center_distance - s.set.radius))
        pos[:, s.num_flap_D]    ~ pos[:, s.num_D] - e_x * flap_length * cos(flap_angle[2]) + e_r_D * flap_length * sin(flap_angle[2]) +
            e_z * (sin(s.α_C) * s.set.radius + (s.set.bridle_center_distance - s.set.radius))
        vel[:, s.num_flap_C]    ~ vel[:, s.num_C] - e_x * flap_length * cos(flap_vel[1]) + e_r_C * flap_length * sin(flap_vel[1])
        vel[:, s.num_flap_D]    ~ vel[:, s.num_D] - e_x * flap_length * cos(flap_vel[2]) + e_r_D * flap_length * sin(flap_vel[2])
        acc[:, s.num_flap_C]    ~ acc[:, s.num_C] - e_x * flap_length * cos(flap_acc[1]) + e_r_C * flap_length * sin(flap_acc[1])
        acc[:, s.num_flap_D]    ~ acc[:, s.num_D] - e_x * flap_length * cos(flap_acc[2]) + e_r_D * flap_length * sin(flap_acc[2])
        segment_length          ~ tether_length  ./ s.set.segments
        mass_tether_particle    ~ mass_per_meter .* segment_length
        damping                 ~ [s.damping / segment_length[1], s.damping / segment_length[2], s.damping*2 / segment_length[3]]
        c_spring                ~ [s.c_spring / segment_length[1], s.c_spring / segment_length[2], s.c_spring*2 / segment_length[3]]
        P_c     ~ 0.5 * (pos[:, s.num_C] + pos[:, s.num_D])
        e_y     ~ (pos[:, s.num_C] - pos[:, s.num_D]) / norm(pos[:, s.num_C] - pos[:, s.num_D])
        e_z     ~ (pos[:, s.num_E] - P_c) / norm(pos[:, s.num_E] - P_c)
        e_x     ~ cross(e_y, e_z)
        e_r_C   ~ (E_c - pos[:, s.num_C]) / norm(E_c - pos[:, s.num_C])
        e_r_D   ~ (E_c - pos[:, s.num_D]) / norm(E_c - pos[:, s.num_D])
        e_te_C  ~ e_x * sin(flap_angle[1]) + e_r_C * cos(flap_angle[1])
        e_te_D  ~ e_x * sin(flap_angle[2]) + e_r_D * cos(flap_angle[2])
        E_c     ~ pos[:, s.num_E] + e_z * (-s.set.bridle_center_distance + s.set.radius) # E_c is the center of the circle shape of the front view of the kite
        rho_kite        ~ calc_rho(s.am, pos[3,s.num_A])
        damping_coeff   ~ max(1.0 - t, 0.0) * s.damping_coeff
    ]

    @variables begin
        winch_force(t)[1:3] # normalized winch forces
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
        distance(t)[1:s.num_A]
        distance_acc(t)
        kite_acc(t)[1:3]
        x_acc(t)
        y_acc(t)
        left_diff(t)
        right_diff(t)
    end
    seqs = [
        seqs
        winch_force     ~ [norm(force[:, i]) for i in 1:3]
        # orientation     ~ orient_euler(s; one_point=false) # TODO: this doesn't work. Calc orientation as observable from camera, with only symbolic vars.
        # upwind_dir      ~ upwind_dir(v_wind_gnd)
        # heading         ~ calc_heading(orientation, elevation, azimuth; upwind_dir)
        heading_y       ~ calc_heading_y(e_x)
        power_angle         ~ (flap_angle[1] + flap_angle[2]) / 2
        power_vel           ~ (flap_vel[1] + flap_vel[2]) / 2
        steering_angle      ~ flap_angle[2] - flap_angle[1]
        steering_vel        ~ flap_vel[2] - flap_vel[1]
        tether_diff         ~ tether_length[2] - tether_length[1]
        tether_diff_vel     ~ tether_vel[2] - tether_vel[1]
        set_diff            ~ set_values[2] - set_values[1]
        # D(set_values)       ~ [0, 0, 0]

        distance_acc        ~ (acc[:, s.num_C] + acc[:, s.num_D] + acc[:, s.num_A]) ⋅ normalize(P_c - pos[:, 3])
        elevation           ~ atan(pos[3, end] / pos[1, end])
        # elevation_vel     ~ vel in direction up
        azimuth             ~ -atan(pos[2, end] / pos[1, end])
        kite_acc            ~ acc[:, s.num_C] + acc[:, s.num_D] + acc[:, s.num_A]
        x_acc               ~ (acc[:, s.num_C] + acc[:, s.num_D] + acc[:, s.num_A]) ⋅ e_x
        y_acc               ~ (acc[:, s.num_C] + acc[:, s.num_D] + acc[:, s.num_A]) ⋅ e_y
        left_diff           ~ tether_length[1] - tether_length[3]
        right_diff          ~ tether_length[2] - tether_length[3]
    ]
    for i in 1:s.num_A
        seqs = [seqs; distance[i] ~ norm(pos[:, i])]
    end
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
        pos(t)[1:3, 1:s.num_A] # left right middle
        vel(t)[1:3, 1:s.num_A] 
        acc(t)[1:3, 1:s.num_A]
        flap_angle(t)[1:2]  # angle left right / C D
        flap_vel(t)[1:2]    # angular vel
        flap_acc(t)[1:2]                # angular acc
        tether_length(t)[1:3]
        tether_vel(t)[1:3]
        tether_acc(t)[1:3]
        segment_length(t)[1:3]
        mass_tether_particle(t)[1:3]
        damping(t)[1:3]
        damping_coeff(t)
        c_spring(t)[1:3]
        P_c(t)[1:3]
        E_c(t)[1:3]
        e_x(t)[1:3]
        e_y(t)[1:3]
        e_z(t)[1:3]
        e_r_C(t)[1:3]
        e_r_D(t)[1:3]
        e_te_C(t)[1:3]
        e_te_D(t)[1:3]
        force(t)[1:3, 1:s.num_A]
        rho_kite(t)
    end
    # Collect the arrays into variables
    pos = collect(pos)
    vel = collect(vel)
    acc = collect(acc)

    deqs = []
    mass_per_meter = s.set.rho_tether * π * (s.set.d_tether/2000.0)^2
    flap_length = s.kite_length_C/4
    @assert flap_length != 0

    [deqs = vcat(deqs, pos[:, i] .~ 0.0) for i in 1:3]
    [deqs = vcat(deqs, D.(pos[:, i]) .~ vel[:, i]) for i in 4:s.num_flap_C-1]
    deqs = [deqs; D(flap_angle)   ~ flap_vel]
    [deqs = vcat(deqs, D.(pos[:, i]) .~ vel[:, i]) for i in s.num_E:s.num_A]
    [deqs = vcat(deqs, vel[:, i] .~ 0.0) for i in 1:3]
    [deqs = vcat(deqs, D.(vel[:, i]) .~ acc[:, i]) for i in 4:s.num_flap_C-1]
    deqs = [deqs; D(flap_vel)   ~ flap_acc]
    [deqs = vcat(deqs, D.(vel[:, i]) .~ acc[:, i]) for i in s.num_E:s.num_A]
    deqs = vcat(deqs, D.(tether_length) .~ tether_vel)
    deqs = vcat(deqs, D.(tether_vel) .~ tether_acc)

    # Compute the masses and forces
    force_eqs = SizedArray{Tuple{3, s.num_A}, Symbolics.Equation}(undef)
    force_eqs[:, :] .= (force[:, :] .~ 0)
    
    seqs = []
    seqs            = scalar_eqs(s, seqs, pos, vel, acc, flap_angle, flap_vel, flap_acc, segment_length, mass_tether_particle, damping, c_spring, 
                        P_c, e_y, e_z, e_x, e_r_C, e_r_D, e_te_C, e_te_D, E_c, rho_kite, damping_coeff, v_wind_gnd, tether_length, tether_vel,
                        mass_per_meter, force, set_values)
    seqs, force_eqs = calc_aero_forces!(s, seqs, force_eqs, force, pos, vel, t, e_x, e_y, e_z, E_c, rho_kite, v_wind, flap_angle)
    seqs, force_eqs = inner_loop!(s, seqs, force_eqs, t, force, pos, vel, segment_length, c_spring, damping, v_wind_gnd)
    
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
    for i in 4:s.num_flap_C-1
        seqs = vcat(seqs, vcat(force_eqs[:, i]))
        seqs = vcat(seqs, acc[:, i] .~ [0.0; 0.0; -G_EARTH] .+ (force[:, i] ./ mass_tether_particle[(i-1)%3+1]) .- damping_coeff * vel[:, i])
    end

    # torque = I * flap_acc
    # flap_acc = torque / (1/3 * (kite_mass/8) * kite_length_c^2)
    # torque = force[:, i] * kite_length_c
    # flap_acc = force[:, i] * kite_length_c / (1/3 * (kite_mass/8) * kite_length_c^2)

    # 1. add all flap + spring + drag forces to flap_C point
    # 2. remove forces not in e_flap_c direction
    # 3. substract forces on point flap_C from point C
    # 4. calculate acceleration from force flap c in e_flap_c direction
    [force_eqs[j, s.num_C] = force[j, s.num_C] ~ force_eqs[j, s.num_C].rhs - force_eqs[j, s.num_flap_C].rhs for j in 1:3]
    [force_eqs[j, s.num_D] = force[j, s.num_D] ~ force_eqs[j, s.num_D].rhs - force_eqs[j, s.num_flap_D].rhs for j in 1:3]
    seqs = [
        seqs
        vcat(force_eqs[:, s.num_flap_C])
        vcat(force_eqs[:, s.num_flap_D])
        flap_acc[1] ~ ((force[:, s.num_flap_C]) ⋅ e_te_C - s.damping * s.flap_damping * flap_vel[1]) * # TODO: add turning drag instead of damping
                    flap_length / (1/3 * (s.set.mass/8) * flap_length^2) - (damping_coeff*200) * flap_vel[1]
        flap_acc[2] ~ ((force[:, s.num_flap_D]) ⋅ e_te_D - s.damping * s.flap_damping * flap_vel[2]) * 
                    flap_length / (1/3 * (s.set.mass/8) * flap_length^2) - (damping_coeff*200) * flap_vel[2]
    ]

    for i in s.num_E:s.num_A
        seqs = vcat(seqs, vcat(force_eqs[:, i]))
        seqs = vcat(seqs, acc[:, i] .~ [0.0; 0.0; -G_EARTH] .+ (force[:, i] ./ s.masses[i]) .- (damping_coeff) * vel[:, i])
    end

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

    if (real) normal_pos_idxs = vcat(4:s.num_flap_C-1, s.num_E)
    else normal_pos_idxs = vcat(4:s.num_flap_C-1, s.num_E, s.num_A) end
    if real
        u0map = [
            # sys.distance[s.num_A] => norm(s.pos[s.num_A])
            [sys.pos[j, s.num_A] => s.pos[s.num_A][j] for j in 1:3]
            [sys.vel[j, i] => norm(s.pos[i]) / norm(s.pos[s.num_A]) * sys.vel[j, s.num_A] 
                for j in 1:3 for i in vcat(4:s.num_flap_C-1, s.num_flap_D+1:s.num_A-1)]
            # [sys.vel[j, i] => norm(s.pos[i]) / norm(s.pos[s.num_A]) * sys.vel[j, s.num_A] 
            #     for j in 1:3 for i in vcat(4:s.num_flap_C-1, s.num_flap_D+1:s.num_A-1)]
            [sys.acc[:, i] => 0 for i in vcat(4:s.num_flap_C-1, s.num_flap_D+1:s.num_A)]

            [sys.flap_vel[j] => 0 for j in 1:2]
            [sys.flap_acc[j] => 0 for j in 1:2]

            sys.gust_factor => 1.0

            # [sys.tether_length[j] => s.tether_lengths[j] for j in 1:3]
            # [sys.tether_vel[j] => [0.01, 0.01, -10][j] for j in 1:3]
            [sys.winch_force[j] => [3.49, 3.49, 57.52][j] for j in 1:3]
            [sys.tether_acc[j] => 0 for j in 1:3]
        ]
    else
        u0map = [
            # sys.distance[s.num_A] => norm(s.pos[s.num_A])
            [sys.pos[j, s.num_A] => s.pos[s.num_A][j] for j in 1:3]
            [sys.vel[j, i] => norm(s.pos[i]) / norm(s.pos[s.num_A]) * sys.vel[j, s.num_A] 
                for j in 1:3 for i in vcat(4:s.num_flap_C-1, s.num_flap_D+1:s.num_A-1)]
            [sys.acc[:, i] => 0 for i in vcat(4:s.num_flap_C-1, s.num_flap_D+1:s.num_A)]

            [sys.flap_vel[j] => 0 for j in 1:2]
            [sys.flap_acc[j] => 0 for j in 1:2]

            sys.gust_factor => 1.0

            [sys.tether_length[j] => s.tether_lengths[j] for j in 1:3]
            # [sys.tether_vel[j] => [1.46, 1.46, 1.46][j] for j in 1:3]
            [sys.tether_acc[j] => 0 for j in 1:3]
        ]
    end
    guesses = [
        [sys.pos[j, i] => s.pos[i][j] for j in 1:3 for i in 1:s.num_A]
        [sys.vel[j, i] => 0 for j in 1:3 for i in 1:s.num_A]

        [sys.flap_angle[j] => 0 for j in 1:2]
        [sys.flap_vel[j] => 0 for j in 1:2]

        [sys.tether_length[j] => s.tether_lengths[j] for j in 1:3]
        [sys.tether_vel[j] => 0 for j in 1:3]

        sys.gust_factor => 1.0
        [sys.kite_acc[j] => 0 for j in 1:3]
        [sys.set_values[j] => s.measure.winch_torque[j] for j in 1:3]
        sys.set_diff => s.measure.winch_torque[2] - s.measure.winch_torque[1]
    ]
    @time prob = ModelingToolkit.InitializationProblem(sys, 0.0, u0map; guesses, fully_determined=true)
    tol = 1e-3
    # @time remake(prob; u0=u0map)
    # @time sol = solve(prob, RobustMultiNewton(); maxiters=10_000, abstol=tol, reltol=tol)
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