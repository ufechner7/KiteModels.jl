
# ==================== mtk model functions ================================================
# Implementation of the three-line model using ModellingToolkit.jl


function calc_acc_speed(motor::AsyncMachine, tether_vel, norm_, set_speed)
    calc_acceleration(motor, tether_vel, norm_; set_speed, set_torque=nothing, use_brake=false) # TODO: add brake setting
end
@register_symbolic calc_acc_speed(motor::AsyncMachine, tether_vel, norm_, set_speed)

function WinchModels.calc_acceleration(wm::TorqueControlledMachine, speed, force; set_torque=nothing, set_speed=nothing, use_brake = false)
    omega      = wm.set.gear_ratio/wm.set.drum_radius * speed
    τ = WinchModels.calc_coulomb_friction(wm) * WinchModels.smooth_sign(omega) + WinchModels.calc_viscous_friction(wm, omega)
    # calculate tau based on the set_torque
    K = 1.0
    tau = set_torque * K
    # calculate tau_total based on the friction
    tau_total = tau + wm.set.drum_radius / wm.set.gear_ratio * force  - τ
    # calculate omega_dot_m based on tau_total and the inertia
    omega_dot_m = tau_total/wm.set.inertia_total
    return wm.set.drum_radius/wm.set.gear_ratio * omega_dot_m
end
function calc_acc_torque(motor::TorqueControlledMachine, tether_vel, norm_, set_torque)
    calc_acceleration(motor, tether_vel, norm_; set_speed=nothing, set_torque, use_brake=false)
end
@register_symbolic calc_acc_torque(motor::TorqueControlledMachine, tether_vel, norm_, set_torque)

function sym_interp(interp::Function, aoa, flap_angle)
    return interp(rad2deg(aoa), rad2deg(flap_angle-aoa))
end
@register_symbolic sym_interp(interp::Function, aoa, flap_angle)


"""
    calc_aero_forces!(s::KPS4_3L, eqs, force_eqs, force, pos, vel, t, e_x, e_y, e_z, rho)

Calculates the aerodynamic forces acting on the kite particles.

Parameters:
- pos:              vector of the particle positions
- vel:              vector of the particle velocities
- rho:              air density [kg/m^3]

Updates the vector s.forces of the first parameter.
"""
function calc_aero_forces!(s::KPS4_3L, eqs, force_eqs, force, pos, vel, t, e_x, e_y, e_z, E_C, rho, v_wind, flap_angle)
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

    eqs = [
        eqs
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
        eqs = [
            eqs
            F[:, i]          ~ E_C + e_y * cos(α) * s.set.radius - e_z * sin(α) * s.set.radius
            e_r[:, i]        ~ (E_C - F[:, i]) / norm(E_C - F[:, i])
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

    
    eqs = [
        eqs
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
    return eqs, force_eqs
end

""" 
    calc_particle_forces!(s::KPS4_3L, eqs, force_eqs, force, pos1, pos2, vel1, vel2, length, c_spring, damping, 
                          rho, i, l_0, k, c, segment, rel_vel, av_vel, norm1, unit_vector, k1, k2, c1, spring_vel,
                          spring_force, v_apparent, v_wind_tether, area, v_app_perp, half_drag_force)

Calculate the drag force and spring force of the tether segment, defined by the parameters pos1, pos2, vel1 and vel2
and distribute it equally on the two particles, that are attached to the segment.
The result is stored in the array s.forces. 
"""
function calc_particle_forces!(s::KPS4_3L, eqs, force_eqs, force, pos1, pos2, vel1, vel2, length, c_spring, 
    damping, rho, i, l_0, k, c, segment, rel_vel, av_vel, norm1, unit_vector, k1, k2, c1, c2, spring_vel, perp_vel,
            spring_force, v_apparent, v_wind_tether, area, v_app_perp, half_drag_force)
    d_tether = s.set.d_tether/1000.0
    eqs = [
        eqs
        i <= s.set.segments*3 ? l_0 ~ length[(i-1) % 3 + 1] : l_0 ~ s.springs[i].length # Unstressed length
        i <= s.set.segments*3 ? k   ~ c_spring[(i-1) % 3 + 1] :
                                k   ~ s.springs[i].c_spring        # Spring constant
        i <= s.set.segments*3 ? c   ~ damping[(i-1) % 3 + 1] : c ~ s.springs[i].damping # Damping coefficient    
        segment      ~ pos1 - pos2 # TODO: all segments have same length and tension
        rel_vel      ~ vel1 - vel2
        av_vel       ~ 0.5 * (vel1 + vel2)
        norm1        ~ norm(segment)
        unit_vector  ~ segment / norm1
        k1           ~ 1.0 * k # compression stiffness kite segments
        k2           ~ 0.1 * k  # compression stiffness tether segments
        c1           ~ 6.0 * c  # damping kite segments
        c2           ~ 0.05 * c  # damping perpendicular
        spring_vel   ~ rel_vel ⋅ unit_vector
        perp_vel     ~ rel_vel .- spring_vel * unit_vector
    ]

    if i >= Base.length(s.springs) - KITE_SPRINGS_3L + 1  # kite springs
        for j in 1:3
            eqs = [
                eqs
                spring_force[j] ~ 
                    (k  * (l_0 - norm1) - c1 * spring_vel) * unit_vector[j] * (1 + smooth_sign_ϵ(norm1 - l_0; s.ϵ)) / 2 +
                    (k1 * (l_0 - norm1) -  c * spring_vel) * unit_vector[j] * (1 - smooth_sign_ϵ(norm1 - l_0; s.ϵ)) / 2
            ]
        end
    else
        for j in 1:3
            eqs = [
                eqs
                spring_force[j] ~
                    ((k  * (l_0 - norm1) - c * spring_vel) * unit_vector[j] - c2 * perp_vel[j]) * (1 + smooth_sign_ϵ(norm1 - l_0; s.ϵ)) / 2 +
                    ((k2 * (l_0 - norm1) - c * spring_vel) * unit_vector[j] - c2 * perp_vel[j]) * (1 - smooth_sign_ϵ(norm1 - l_0; s.ϵ)) / 2
            ]
        end
    end
    eqs = [
        eqs
        v_apparent       ~ v_wind_tether - av_vel
        i >= s.num_flap_C ?
            area             ~ norm1 * d_tether * 10 : # 10 is the number of parallel lines in the bridle system
            area             ~ norm1 * d_tether * (1 + (i%3 == 0)) # double area for middle tether
        v_app_perp       ~ v_apparent - (v_apparent ⋅ unit_vector) * unit_vector
        half_drag_force  ~ (0.25 * rho * s.set.cd_tether * norm(v_app_perp) * area) .* v_app_perp
    ]

    for j in 1:3
        force_eqs[j, s.springs[i].p1] = 
            (force[j, s.springs[i].p1] ~ force_eqs[j, s.springs[i].p1].rhs + (half_drag_force[j] + spring_force[j]))
        force_eqs[j, s.springs[i].p2] = 
            (force[j, s.springs[i].p2] ~ force_eqs[j, s.springs[i].p2].rhs + (half_drag_force[j] - spring_force[j]))
    end
    
    return eqs, force_eqs
end



"""
    inner_loop!(s::KPS4_3L, eqs, force_eqs, t, force, pos, vel, length, c_spring, damping, v_wind_gnd)

Calculate the forces, acting on all particles.

Output:length
- s.forces
- s.v_wind_tether
"""
@inline function inner_loop!(s::KPS4_3L, eqs, force_eqs, t, force, pos, vel, length, c_spring, damping, v_wind_gnd)
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
    end
    
    for i in eachindex(s.springs)
        p1 = Int(s.springs[i].p1)  # First point nr.
        p2 = Int(s.springs[i].p2)
        eqs = [
            eqs
            height[i]           ~ max(0.0, 0.5 * (pos[:, p1][3] + pos[:, p2][3]))
            rho[i]              ~ calc_rho(s.am, height[i])
            v_wind_tether[:, i] ~ AtmosphericModels.calc_wind_factor(s.am, height[i], s.set.profile_law) * v_wind_gnd
        ]

        eqs, force_eqs = calc_particle_forces!(s, eqs, force_eqs, force, pos[:, p1], pos[:, p2], vel[:, p1], 
                          vel[:, p2], length, c_spring, damping, rho[i], i, l_0[i], k[i], c[i], segment[:, i], 
                          rel_vel[:, i], av_vel[:, i], norm1[i], unit_vector[:, i], k1[i], k2[i], c1[i], c2[i], spring_vel[i],
                          perp_vel[:, i], spring_force[:, i], v_apparent[:, i], v_wind_tether[:, i], area[i], v_app_perp[:, i], 
                          half_drag_force[:, i])
    end

    return eqs, force_eqs
end

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

function model!(s::KPS4_3L, pos, vel, e_x, e_y, e_z, flap_angle, eqs=[])
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
        set_values(t)[1:3] = zeros(3) # left right middle
        acc(t)[1:3, 1:s.num_A]
        flap_vel(t)[1:2]     = zeros(2) # angular vel
        flap_acc(t)[1:2]                # angular acc
        tether_length(t)[1:3]  = s.tether_lengths
        tether_vel(t)[1:3] = zeros(3)
        tether_acc(t)[1:3]
        segment_length(t)[1:3]
        mass_tether_particle(t)[1:3]
        damping(t)[1:3]
        c_spring(t)[1:3]
        P_c(t)[1:3]
        E_C(t)[1:3]
        e_r_C(t)[1:3]
        e_r_D(t)[1:3]
        e_te_C(t)[1:3]
        e_te_D(t)[1:3]
        force(t)[1:3, 1:s.num_A]
        rho_kite(t)
        winch_force(t)[1:3] # normalized winch forces
        heading(t)
        heading_y(t)
        power_angle(t) # average flap angle
        power_vel(t)
        steering_angle(t) # difference between left and right flap angle
        steering_vel(t)
        tether_diff(t)
        tether_diff_vel(t)
        set_diff(t)
    end
    # Collect the arrays into variables
    # acc = collect(acc)

    mass_per_meter = s.set.rho_tether * π * (s.set.d_tether/2000.0)^2

    # Compute the masses and forces
    force_eqs = SizedArray{Tuple{3, s.num_A}, Symbolics.Equation}(undef)
    force_eqs[:, :] .= (force[:, :] .~ 0)

    flap_length = s.kite_length_C/4
    eqs = [
        eqs
        vel[:, s.num_flap_C]    ~ vel[:, s.num_C] - e_x * flap_length * cos(flap_vel[1]) + e_r_C * flap_length * sin(flap_vel[1])
        vel[:, s.num_flap_D]    ~ vel[:, s.num_D] - e_x * flap_length * cos(flap_vel[2]) + e_r_D * flap_length * sin(flap_vel[2])
        acc[:, s.num_flap_C]    ~ acc[:, s.num_C] - e_x * flap_length * cos(flap_acc[1]) + e_r_C * flap_length * sin(flap_acc[1])
        acc[:, s.num_flap_D]    ~ acc[:, s.num_D] - e_x * flap_length * cos(flap_acc[2]) + e_r_D * flap_length * sin(flap_acc[2])
        segment_length          ~ tether_length  ./ s.set.segments
        mass_tether_particle    ~ mass_per_meter .* segment_length
        damping                 ~ [s.damping / segment_length[1], s.damping / segment_length[2], s.damping*2 / segment_length[3]]
        c_spring                ~ [s.c_spring / segment_length[1], s.c_spring / segment_length[2], s.c_spring*2 / segment_length[3]]
        P_c     ~ 0.5 * (pos[:, s.num_C] + pos[:, s.num_D])
        e_r_C   ~ (E_C - pos[:, s.num_C]) / norm(E_C - pos[:, s.num_C])
        e_r_D   ~ (E_C - pos[:, s.num_D]) / norm(E_C - pos[:, s.num_D])
        e_te_C  ~ e_x * sin(flap_angle[1]) + e_r_C * cos(flap_angle[1])
        e_te_D  ~ e_x * sin(flap_angle[2]) + e_r_D * cos(flap_angle[2])
        E_C     ~ pos[:, s.num_E] + e_z * (-s.set.bridle_center_distance + s.set.radius) # E_C is the center of the circle shape of the front view of the kite
        rho_kite        ~ calc_rho(s.am, pos[:, s.num_A][3])
        winch_force     ~ [norm(force[:, i]) for i in 1:3]
        heading         ~ calc_heading(e_x, pos[:, s.num_E])
        heading_y       ~ calc_heading_y(e_x)
        power_angle         ~ (flap_angle[1] + flap_angle[2]) / 2
        power_vel           ~ (flap_vel[1] + flap_vel[2]) / 2
        steering_angle      ~ flap_angle[2] - flap_angle[1]
        steering_vel        ~ flap_vel[2] - flap_vel[1]
        tether_diff         ~ tether_length[2] - tether_length[1]
        tether_diff_vel     ~ tether_vel[2] - tether_vel[1]
        set_diff            ~ set_values[2] - set_values[1]
    ]

    if s.torque_control
        eqs = [eqs; tether_acc ~ [calc_acc_torque(s.motors[i], tether_vel[i], norm(force[:, (i-1)%3+1]),
            set_values[i]) for i in 1:3]]
    else
        eqs = [eqs; tether_acc ~ [calc_acc_speed(s.motors[i], tether_vel[i], norm(force[:, (i-1)%3+1]),
            set_values[i]) for i in 1:3]]
    end
    for i in 1:3
        eqs = [eqs; vcat(force_eqs[:, i])]
        eqs = [eqs; acc[:, i] ~ zeros(3)]
    end
    for i in 4:s.num_flap_C-1
        eqs = [eqs; vcat(force_eqs[:, i])]
        eqs = [eqs; acc[:, i] ~ [0.0; 0.0; -G_EARTH] .+ (force[:, i] ./ mass_tether_particle[(i-1)%3+1])]
    end

    eqs, force_eqs = calc_aero_forces!(s, eqs, force_eqs, force, pos, vel, t, e_x, e_y, e_z, E_C, rho_kite, v_wind, flap_angle)
    eqs, force_eqs = inner_loop!(s, eqs, force_eqs, t, force, pos, vel, segment_length, c_spring, damping, v_wind_gnd)
    
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
    eqs = [
        eqs
        vcat(force_eqs[:, s.num_flap_C])
        vcat(force_eqs[:, s.num_flap_D])
        flap_acc[1] ~ ((force[:, s.num_flap_C]) ⋅ e_te_C - s.damping * s.flap_damping * flap_vel[1]) * # TODO: add turning drag instead of damping
                    flap_length / (1/3 * (s.set.mass/8) * flap_length^2)
        flap_acc[2] ~ ((force[:, s.num_flap_D]) ⋅ e_te_D - s.damping * s.flap_damping * flap_vel[2]) * 
                    flap_length / (1/3 * (s.set.mass/8) * flap_length^2)
    ]

    for i in s.num_E:s.num_A
        eqs = [eqs; vcat(force_eqs[:, i])]
        eqs = [eqs; acc[:, i] ~ [0.0; 0.0; -G_EARTH] .+ (force[:, i] ./ s.masses[i])]
    end
    return eqs, (acc, tether_length, tether_vel, tether_acc, flap_angle, flap_vel, flap_acc, set_values, heading_y)
end



function ode_model!(s::KPS4_3L, pos_, vel_)
    pos_, vel_ = convert_pos_vel(s, pos_, vel_)
    @variables begin
        pos(t)[1:3, 1:s.num_A] = pos_ # left right middle
        vel(t)[1:3, 1:s.num_A] = vel_
        flap_angle(t)[1:2]   = deg2rad(s.set.alpha_zero) # angle left right / C D
        e_x(t)[1:3]
        e_y(t)[1:3]
        e_z(t)[1:3]
    end
    pos = collect(pos)
    vel = collect(vel)

    eqs, (acc, tether_length, tether_vel, tether_acc, flap_angle, flap_vel, flap_acc, set_values, heading_y) = 
        model!(s, pos, vel)

    eqs = [
        eqs
        e_y     ~ (pos[:, s.num_C] - pos[:, s.num_D]) / norm(pos[:, s.num_C] - pos[:, s.num_D])
        e_z     ~ (pos[:, s.num_E] - P_c) / norm(pos[:, s.num_E] - P_c)
        e_x     ~ e_y × e_z
        pos[:, s.num_flap_C]    ~ pos[:, s.num_C] - e_x * flap_length * cos(flap_angle[1]) + e_r_C * flap_length * sin(flap_angle[1])
        pos[:, s.num_flap_D]    ~ pos[:, s.num_D] - e_x * flap_length * cos(flap_angle[2]) + e_r_D * flap_length * sin(flap_angle[2])
    ]

    [eqs = [eqs; pos[:, i] ~ 0.0] for i in 1:3]
    [eqs = [eqs; D.(pos[:, i]) ~ vel[:, i]] for i in 4:s.num_flap_C-1]
    eqs = [eqs; D(flap_angle)   ~ flap_vel]
    [eqs = [eqs; D.(pos[:, i]) ~ vel[:, i]] for i in s.num_E:s.num_A]
    [eqs = [eqs; vel[:, i] ~ 0.0] for i in 1:3]
    [eqs = [eqs; D.(vel[:, i]) ~ acc[:, i]] for i in 4:s.num_flap_C-1]
    eqs = [eqs; D(flap_vel)   ~ flap_acc]
    [eqs = [eqs; D.(vel[:, i]) ~ acc[:, i]] for i in s.num_E:s.num_A]
    eqs = [eqs; D.(tether_length) ~ tether_vel]
    eqs = [eqs; D.(tether_vel) ~ tether_acc]

    @variables turn_rate_y(t)
    eqs = [eqs; turn_rate_y     ~ D(heading_y)]

    @named sys = ODESystem(Symbolics.scalarize.(reduce(vcat, Symbolics.scalarize.(eqs))), t)
    return sys, collect(set_values)
end


function get_kite_particles(s::KPS4_3L, eqs, pos, e_x, e_y, e_z, flap_angle)
    width, radius, middle_length, tip_length, bridle_center_distance = 
        s.set.width, s.set.radius, s.set.middle_length, s.set.tip_length, s.set.bridle_center_distance
    flap_length = s.kite_length_C/4
    @variables begin
        α_0
        α_C
        α_D
        kite_length_C
        E_c[1:3]
        P_c[1:3]
        e_r_C[1:3]
        e_r_D[1:3]
    end
    eqs = [
        eqs
        α_0 ~ pi/2 - width/2/radius
        α_C ~ α_0 + width*(-2*tip_length + sqrt(2*middle_length^2 + 2*tip_length^2)) /
            (4*(middle_length - tip_length)) / radius
        α_D ~ π - α_C
    
        E_c             ~ pos[:, s.num_E] + e_z * (-bridle_center_distance + radius) # E at center of circle on which the kite shape lies
        pos[:, s.num_C] ~ E_c + e_y*cos(α_C)*radius - e_z*sin(α_C)*radius
        pos[:, s.num_D] ~ E_c + e_y*cos(α_D)*radius - e_z*sin(α_D)*radius
    
        kite_length_C   ~ tip_length + (middle_length-tip_length) * (α_C - α_0) / (π/2 - α_0)
        P_c             ~ (pos[:, s.num_C]+pos[:, s.num_D])./2
        pos[:, s.num_A] ~ P_c - e_x*(kite_length_C*(3/4 - 1/4))

        e_r_C ~ (E_c - pos[:, s.num_C]) / norm(E_c - pos[:, s.num_C])
        e_r_D ~ (E_c - pos[:, s.num_D]) / norm(E_c - pos[:, s.num_D])    
        pos[:, s.num_flap_C] ~ pos[:, s.num_C] - e_x * flap_length * cos(flap_angle[1]) + e_r_C * flap_length * sin(flap_angle[1])
        pos[:, s.num_flap_D] ~ pos[:, s.num_D] - e_x * flap_length * cos(flap_angle[2]) + e_r_D * flap_length * sin(flap_angle[2])
    ]
    return eqs
end


"""
measurements:
    - middle tether length, unstretched
    - azimuth
    - elevation
    - heading
    - first and second order derivatives

    no need for sensors on the kite, just one camera pointing at the kite
    maybe need wind sensor, but middle tether combined with ground wind speed is basically a wind sensor

unknowns:
    - left and right tether length
    - wind speed at kite position
    - tether bend
    - flap angle
    - stretched middle tether length / kite distance
    - kite to middle tether angle

building kite points:
    - set point E at kite_distance
    - build middle tether with bend
    - find e_x, e_y, e_z using middle tether bend
    - build rest of kite points using these axes and using flap_angles
    - build left and right tether with left and right bend
    - 
"""
function nonlin_model!(s::KPS4_3L, pos_, vel_)
    # measurements
    @parameters begin
        m_azimuth = 0.0
        m_elevation = 85
        m_tether_length = 50.0
        m_tether_vel[1:3] = zeros(3)
        m_tether_acc[1:3] = zeros(3)
        m_heading_y = 0.0
        m_turn_rate_y = 0.0
        m_turn_acc_y = 0.0
    end
    @variables begin
        pos[1:3, 1:s.num_A]     = ones(3, s.num_A) # left right middle
        vel[1:3, 1:s.num_A]     = zeros(3, s.num_A)
        r[1:3]                  = ones(3)
        kite_distance[1:3]      = s.tether_lengths  # distance from ground station to kite tether connection point
        flap_angle[1:2]         = zeros(2)
        tether_rot[1:3]         = zeros(3)
        t_y[1:3, 1:3]           = repeat([0, 1, 0]', 3, 1) # tether y vector in ENU frame
        t_x[1:3, 1:3]           = repeat([1, 0, 0]', 3, 1) 
        t_z[1:3, 1:3]           = repeat([0, 0, 1]', 3, 1) 
        e_x[1:3]                = repeat([0, 1, 0]', 3, 1) 
        e_y[1:3]                = repeat([1, 0, 0]', 3, 1) 
        e_z[1:3]                = repeat([0, 0, 1]', 3, 1) 
        ls_y[1:3]               = repeat([0, 1, 0]', 3, 1)  # last tether segment y vector in NED frame
        ls_x[1:3]               = repeat([-1, 0, 0]', 3, 1) 
        ls_z[1:3]               = repeat([0, 0, -1]', 3, 1) 
        kite_angle[1:2]         = zeros(2) # angles between last tether segment and kite reference frame
    end

    # calculate segment length
    # segment length is known because of the measured tether force
    # with the known segment length and wind speed there should be 1 solution for tether bending
    eqs = [
        t_z[:, 3] ~ rotate_in_yx((rotate_in_xz([0.0, 1.0, 0.0], m_elevation)), m_azimuth)
        # build kite points from bend + heading + kite_distance
        pos[:, s.num_E] ~ t_z[:, 3] * kite_distance[3]
        ls_z ~ (pos[:, s.num_E-3] - pos[:, s.num_E]) / norm(pos[:, s.num_E-3] - pos[:, s.num_E])
        calc_heading_y(ls_x) ~ m_heading_y # TODO: check if this works
        ls_x ~ ls_y × ls_z
        ls_y ~ ls_z × ls_x
        e_x ~ rotate_v_k_angle(rotate_v_k_angle(ls_x, ls_y, kite_angle[1]), ls_x, kite_angle[2])
        e_z ~ rotate_v_k_angle(rotate_v_k_angle(ls_z, ls_y, kite_angle[1]), ls_x, kite_angle[2])
        e_y ~ rotate_v_k_angle(rotate_v_k_angle(ls_y, ls_y, kite_angle[1]), ls_x, kite_angle[2])
    ]

    # generate kite points
    eqs = get_kite_particles(s, eqs, pos, e_x, e_y, e_z, flap_angle)

    eqs = [
        eqs
        t_z[:, 1] ~ pos[:, s.num_flap_C] / norm(pos[:, s.num_flap_C])
        t_z[:, 2] ~ pos[:, s.num_flap_D] / norm(pos[:, s.num_flap_D])
        kite_distance[1] ~ norm(pos[:, s.num_flap_C])
        kite_distance[2] ~ norm(pos[:, s.num_flap_D])
    ]    
    for i in 1:3
        eqs = [
            eqs
            t_y[:, i] ~ e_x × t_z[:, i]
            t_x[:, i] ~ t_y[:, i] × t_z[:, i]    
        ]
    end

    # draw tethers
    @variables begin
        α[1:3]
        γ[1:3, 1:s.set.segments-1]
        local_z[1:3, 1:s.set.segments-1]
        local_x[1:3, 1:s.set.segments-1]
        local_y[1:3, 1:s.set.segments-1]
        h[1:3]
        r[1:3]
    end
    for i in 1:3
        eqs = [
            eqs
            h[i] ~ kite_distance[i]/(2tan(α[i]))
            r[i] ~ kite_distance[i]/(2sin(α[i]))
        ]
    end
    for (i, j) in enumerate(range(4, step=3, length=s.set.segments-1))
        eqs = [
            eqs
            γ[:, i]         ~ -α + 2α*i / s.set.segments
            local_z[:, i]   ~ (kite_distance/2 .+ r .* sin.(γ[:, i]))
            local_x[:, i]   ~ (r .* cos.(γ[:, i]) - h) .* sin.(tether_rot)
            local_y[:, i]   ~ (r .* cos.(γ[:, i]) - h) .* cos.(tether_rot)
            ]
        for k in 1:3
            eqs = [
                eqs
                pos[:, j+k-1] ~ local_z[k, i] * t_z[:, k] + local_x[k, i] * t_x[:, k] + local_y[k, i] * t_y[:, k]
            ]
        end
    end

    eqs, (acc, tether_length, tether_vel, tether_acc, flap_vel, flap_acc, set_values, heading_y) = 
        model!(s, pos, vel, e_x, e_y, e_z, flap_angle, eqs)

    
    [eqs =  [eqs; zeros(3) ~ pos[:, i]] for i in 1:3]
    [eqs =  [eqs; zeros(3) ~ vel[:, i]] for i in 1:3]
    [eqs =  [eqs; zeros(3) ~ vel[:, i]] for i in 4:s.num_flap_C-1]
    [eqs =  [eqs; zeros(3) ~ vel[:, i]] for i in s.num_E:s.num_A]
    eqs =   [eqs; zeros(2) ~ flap_vel]
    # TODO: spring forces in kite segments ~ 0.0
    [eqs =  [eqs; zeros(3) ~ acc[:, i]] for i in 1:3]
    [eqs =  [eqs; zeros(3) ~ acc[:, i]] for i in 4:s.num_flap_C-1]
    [eqs =  [eqs; zeros(3) ~ acc[:, i]] for i in s.num_E:s.num_A]
    eqs =   [eqs; zeros(2) ~ flap_acc]
    # eqs =   [eqs; 0 ~ tether_vel]
    # eqs =   [eqs; 0 ~ tether_acc]

    eqs = [
        eqs
        tether_length[3] ~ m_tether_length
        tether_acc ~ m_tether_acc
        tether_vel ~ m_tether_vel
    ]
    # return eqs

    eqs = Symbolics.scalarize.(reduce(vcat, Symbolics.scalarize.(eqs)))
    @mtkbuild nsys = NonlinearSystem(eqs) fully_determined = false  # TODO: optimalizationproblem
    return nsys
end