# Implementation of the three-line model using ModellingToolkit.jl

function calc_acc_speed(tether_speed::SimFloat, norm_::SimFloat, set_speed::SimFloat)
    calc_acceleration(AsyncMachine(se("system_3l.yaml")), tether_speed, norm_; set_speed, set_torque=nothing, use_brake=false)
end
@register_symbolic calc_acc_speed(tether_speed, norm_, set_speed)

function calc_acc_torque(tether_speed::SimFloat, norm_::SimFloat, set_torque::SimFloat)
    calc_acceleration(TorqueControlledMachine(se("system_3l.yaml")), tether_speed, norm_; set_speed=nothing, set_torque, use_brake=false)
end
@register_symbolic calc_acc_torque(tether_speed, norm_, set_torque)

"""
    calc_aero_forces!(s::KPS4_3L, eqs2, force_eqs, force, pos, vel, t, e_x, e_y, e_z, rho)

Calculates the aerodynamic forces acting on the kite particles.

Parameters:
- pos:              vector of the particle positions
- vel:              vector of the particle velocities
- rho:              air density [kg/m^3]

Updates the vector s.forces of the first parameter.
"""
function calc_aero_forces_mtk!(s::KPS4_3L, eqs2, force_eqs, force, pos, vel, t, e_x, e_y, e_z, rho, v_wind)
    n = s.set.aero_surfaces
    @variables begin
        δ_left(t)
        δ_right(t)
        E_c(t)[1:3]
        v_cx(t)[1:3]
        v_dx(t)[1:3]
        v_dy(t)[1:3]
        v_dz(t)[1:3]
        v_cy(t)[1:3]
        v_cz(t)[1:3]
        y_lc(t)
        y_ld(t)
    end

    eqs2 = [
        eqs2
        δ_left  ~ (pos[:, s.num_E-2] - pos[:,s.num_C]) ⋅ e_z
        δ_right ~ (pos[:, s.num_E-1] - pos[:,s.num_D]) ⋅ e_z
        # in the aero calculations, E_c is the center of the circle shape on which the kite lies
        E_c     ~ pos[:, s.num_E] + e_z * (-s.set.bridle_center_distance + s.set.radius)
        v_cx    ~ (vel[:, s.num_C] ⋅ e_x) * e_x
        v_dx    ~ (vel[:, s.num_D] ⋅ e_x) * e_x
        v_dy    ~ (vel[:, s.num_D] ⋅ e_y) * e_y
        v_dz    ~ (vel[:, s.num_D] ⋅ e_z) * e_z
        v_cy    ~ (vel[:, s.num_C] ⋅ e_y) * e_y
        v_cz    ~ (vel[:, s.num_C] ⋅ e_z) * e_z
        y_lc    ~  norm(pos[:, s.num_C] - 0.5 * (pos[:, s.num_C] + pos[:, s.num_D]))
        y_ld    ~ -norm(pos[:, s.num_D] - 0.5 * (pos[:, s.num_C] + pos[:, s.num_D]))
    ]

    # integrating loop variables
    @variables begin
        F(t)[1:3, 1:2n]
        e_r(t)[1:3, 1:2n]
        y_l(t)[1:2n]
        v_kite(t)[1:3, 1:2n]
        v_a(t)[1:3, 1:2n]
        e_drift(t)[1:3, 1:2n]
        v_a_xr(t)[1:3, 1:2n]
        aoa(t)[1:n*2]
        dL_dα(t)[1:3, 1:2n]
        dD_dα(t)[1:3, 1:2n]
        L_C(t)[1:3]
        L_D(t)[1:3]
        D_C(t)[1:3]
        D_D(t)[1:3]
        F_steering_c(t)[1:3]
        F_steering_d(t)[1:3]
        d(t)[1:2n]
    end
    l_c_eq = SizedArray{Tuple{3}, Symbolics.Equation}(collect(L_C .~ 0))
    l_d_eq = SizedArray{Tuple{3}, Symbolics.Equation}(collect(L_D .~ 0))
    d_c_eq = SizedArray{Tuple{3}, Symbolics.Equation}(collect(D_C .~ 0))
    d_d_eq = SizedArray{Tuple{3}, Symbolics.Equation}(collect(D_D .~ 0))
    kite_length = zeros(MVector{2n, SimFloat})
    α           = zero(SimFloat)
    α_0         = zero(SimFloat)
    α_middle    = zero(SimFloat)
    dα          = zero(SimFloat)
    # Calculate the lift and drag
    α_0         = π/2 - s.set.width/2/s.set.radius
    α_middle    = π/2
    dα          = (α_middle - α_0) / n
    for i in 1:n*2
        if i <= n
            α = α_0 + -dα/2 + i * dα
        else
            α = pi - (α_0 + -dα/2 + (i-n) * dα)
        end

        eqs2 = [
            eqs2
            F[:, i]          ~ E_c + e_y * cos(α) * s.set.radius - e_z * sin(α) * s.set.radius
            e_r[:, i]        ~ (E_c - F[:, i]) / norm(E_c - F[:, i])
            y_l[i]           ~ cos(α) * s.set.radius
            α < π/2 ?
                v_kite[:, i] ~ ((v_cx - v_dx) / (y_lc - y_ld) * (y_l[i] - y_ld) + v_dx) + v_cy + v_cz :
                v_kite[:, i] ~ ((v_cx - v_dx) / (y_lc - y_ld) * (y_l[i] - y_ld) + v_dx) + v_dy + v_dz
            v_a[:,i]         ~ v_wind .- v_kite[:,i]
            e_drift[:, i]    ~ (e_r[:, i] × e_x)
            v_a_xr[:, i]     ~ v_a[:, i] .- (v_a[:, i] ⋅ e_drift[:, i]) .* e_drift[:, i]
        ]
        if α < π/2
            kite_length[i] = (s.set.tip_length + (s.set.middle_length-s.set.tip_length) * α 
                              * s.set.radius/(0.5*s.set.width))
        else
            kite_length[i] = (s.set.tip_length + (s.set.middle_length-s.set.tip_length) * (π-α) 
                              * s.set.radius/(0.5*s.set.width))
        end
        eqs2 = [
            eqs2
            α < s.α_l ?
                d[i]    ~ δ_left :
            α > s.α_r ?
                d[i]    ~ δ_right :
                d[i]    ~ (δ_right - δ_left) / (s.α_r - s.α_l) * (α - s.α_l) + (δ_left)
            aoa[i]      ~ -asin((v_a_xr[:, i] / norm(v_a_xr[:, i])) ⋅ e_r[:, i]) + 
                           asin(clamp(d[i] / kite_length[i], -1.0, 1.0))
            dL_dα[:, i] ~ 0.5 * rho * (norm(v_a_xr[:, i]))^2 * s.set.radius * kite_length[i] * rad_cl_mtk(aoa[i]) * 
                                ((v_a_xr[:, i] × e_drift[:, i]) / norm(v_a_xr[:, i] × e_drift[:, i]))
            dD_dα[:, i] ~ 0.5 * rho * norm(v_a_xr[:, i]) * s.set.radius * kite_length[i] * rad_cd_mtk(aoa[i]) * 
                                v_a_xr[:,i] # the sideways drag cannot be calculated with the C_d formula
        ]
        if i <= n
            [l_c_eq[j] = (L_C[j] ~ l_c_eq[j].rhs + dL_dα[j, i] * dα) for j in 1:3]
            [d_c_eq[j] = (D_C[j] ~ d_c_eq[j].rhs + dD_dα[j, i] * dα) for j in 1:3]
        else 
            [l_d_eq[j] = (L_D[j] ~ l_d_eq[j].rhs + dL_dα[j, i] * dα) for j in 1:3]
            [d_d_eq[j] = (D_D[j] ~ d_d_eq[j].rhs + dD_dα[j, i] * dα) for j in 1:3]
        end
    end
    
    eqs2 = [
        eqs2
        l_c_eq
        d_c_eq
        l_d_eq
        d_d_eq
        F_steering_c ~ ((0.2 * (L_C ⋅ -e_z)) * -e_z)
        F_steering_d ~ ((0.2 * (L_D ⋅ -e_z)) * -e_z)
    ]
    
    force_eqs[:,s.num_C]   .= (force[:,s.num_C]   .~ (L_C + D_C) - F_steering_c) 
    force_eqs[:,s.num_D]   .= (force[:,s.num_D]   .~ (L_D + D_D) - F_steering_d)
    force_eqs[:,s.num_E-2] .= (force[:,s.num_E-2] .~ F_steering_c)
    force_eqs[:,s.num_E-1] .= (force[:,s.num_E-1] .~ F_steering_d)
    return eqs2, force_eqs
end

""" 
    calc_particle_forces!(s::KPS4_3L, eqs2, force_eqs, force, pos1, pos2, vel1, vel2, length, c_spring, damping, 
                          rho, i, l_0, k, c, segment, rel_vel, av_vel, norm1, unit_vector, k1, k2, c1, spring_vel,
                          spring_force, v_apparent, v_wind_tether, area, v_app_perp, half_drag_force)

Calculate the drag force and spring force of the tether segment, defined by the parameters pos1, pos2, vel1 and vel2
and distribute it equally on the two particles, that are attached to the segment.
The result is stored in the array s.forces. 
"""
function calc_particle_forces_mtk!(s::KPS4_3L, eqs2, force_eqs, force, pos1, pos2, vel1, vel2, length, c_spring, 
    damping, rho, i, l_0, k, c, segment, rel_vel, av_vel, norm1, unit_vector, k1, k2, c1, spring_vel,
            spring_force, v_apparent, v_wind_tether, area, v_app_perp, half_drag_force, stiffness_factor)
    d_tether = s.set.d_tether/1000.0
    eqs2 = [
        eqs2
        i <= s.set.segments*3 ? l_0 ~ length[(i-1) % 3 + 1] : l_0 ~ s.springs[i].length # Unstressed length
        i <= s.set.segments*3 ? k   ~ c_spring[(i-1) % 3 + 1] * stiffness_factor :
                                k   ~ s.springs[i].c_spring * stiffness_factor        # Spring constant
        i <= s.set.segments*3 ? c   ~ damping[(i-1) % 3 + 1] : c ~ s.springs[i].damping # Damping coefficient    
        segment     .~ pos1 - pos2
        rel_vel     .~ vel1 - vel2
        av_vel      .~ 0.5 * (vel1 + vel2)
        norm1        ~ norm(segment)
        unit_vector .~ segment / norm1
        k1           ~ 0.25 * k # compression stiffness kite segments
        k2           ~ 0.1 * k  # compression stiffness tether segments
        c1           ~ 6.0 * c  # damping kite segments
        spring_vel  .~ rel_vel ⋅ unit_vector
    ]

    if i >= s.num_E-2  # kite springs
        for j in 1:3
            eqs2 = [
                eqs2
                spring_force[j] ~ ifelse(
                    (norm1 - l_0) > 0.0,
                    (k  * (l_0 - norm1) - c1 * spring_vel) * unit_vector[j],
                    (k1 * (l_0 - norm1) -  c * spring_vel) * unit_vector[j]
                )
            ]
        end
    else
        for j in 1:3
            eqs2 = [
                eqs2
                spring_force[j] ~ ifelse(
                    (norm1 - l_0) > 0.0,
                    (k  * (l_0 - norm1) - c * spring_vel) * unit_vector[j],
                    (k2 * (l_0 - norm1) - c * spring_vel) * unit_vector[j]
                    )
            ]
        end
    end
    eqs2 = [
        eqs2
        v_apparent       ~ v_wind_tether - av_vel
        area             ~ norm1 * d_tether
        v_app_perp       ~ v_apparent - v_apparent ⋅ unit_vector * unit_vector
        half_drag_force .~ (0.25 * rho * s.set.cd_tether * norm(v_app_perp) * area) .* v_app_perp
    ]

    for j in 1:3
        force_eqs[j, s.springs[i].p1] = 
            (force[j,s.springs[i].p1] ~ force_eqs[j, s.springs[i].p1].rhs + (half_drag_force[j] + spring_force[j]))
        force_eqs[j, s.springs[i].p2] = 
            (force[j,s.springs[i].p2] ~ force_eqs[j, s.springs[i].p2].rhs + (half_drag_force[j] - spring_force[j]))
    end
    
    return eqs2, force_eqs
end



"""
    inner_loop_mtk!(s::KPS4_3L, eqs2, force_eqs, t, force, pos, vel, length, c_spring, damping, v_wind_gnd)

Calculate the forces, acting on all particles.

Output:length
- s.forces
- s.v_wind_tether
"""
@inline function inner_loop_mtk!(s::KPS4_3L, eqs2, force_eqs, t, force, pos, vel, length, c_spring, damping, v_wind_gnd, stiffness_factor)
    @variables begin
        height(t)[eachindex(s.springs)]
        rho(t)[eachindex(s.springs)]
        v_wind_tether(t)[1:3, eachindex(s.springs)]

        l_0(t)[eachindex(s.springs)]
        k(t)[eachindex(s.springs)]
        c(t)[eachindex(s.springs)]
        segment(t)[1:3,eachindex(s.springs)]
        rel_vel(t)[1:3, eachindex(s.springs)]
        av_vel(t)[1:3, eachindex(s.springs)] 
        norm1(t)[eachindex(s.springs)]
        unit_vector(t)[1:3, eachindex(s.springs)]
        k1(t)[eachindex(s.springs)]
        k2(t)[eachindex(s.springs)]
        c1(t)[eachindex(s.springs)]
        spring_vel(t)[eachindex(s.springs)]
        spring_force(t)[1:3, eachindex(s.springs)]
        v_apparent(t)[1:3, eachindex(s.springs)]
        area(t)[eachindex(s.springs)]
        v_app_perp(t)[1:3, eachindex(s.springs)]
        half_drag_force(t)[1:3, eachindex(s.springs)]
    end
    
    for i in eachindex(s.springs)
        p1 = s.springs[i].p1  # First point nr.
        p2 = s.springs[i].p2
        eqs2 = [
            eqs2
            height[i]           ~ 0.5 * (pos[:, p1][3] + pos[:, p2][3])
            rho[i]              ~ calc_rho(s.am, height[i])
            v_wind_tether[:, i] ~ calc_wind_factor(s.am, height[i]) * v_wind_gnd
        ]

        # TODO: @assert height > 0
        eqs2, force_eqs = calc_particle_forces_mtk!(s, eqs2, force_eqs, force, pos[:, p1], pos[:, p2], vel[:, p1], 
                          vel[:, p2], length, c_spring, damping, rho[i], i, l_0[i], k[i], c[i], segment[:, i], 
                          rel_vel[:, i], av_vel[:, i], norm1[i], unit_vector[:, i], k1[i], k2[i], c1[i], spring_vel[i],
                          spring_force[:, i], v_apparent[:,i], v_wind_tether[:, i], area[i], v_app_perp[:, i], 
                          half_drag_force[:, i], stiffness_factor)
    end

    return eqs2, force_eqs
end

function update_pos!(s, integrator)
    pos = s.get_pos(integrator)
    s.steering_pos     .= s.get_steering_pos(integrator)
    [s.pos[i]          .= pos[:,i] for i in 1:s.num_A]
    s.veld[s.num_E-2]  .= s.get_line_acc(integrator)
    s.vel_kite         .= s.get_kite_vel(integrator)
    winch_forces        = s.get_winch_forces(integrator)
    [s.winch_forces[i] .= (winch_forces[:,i]) for i in 1:3]
    s.tether_lengths   .= s.get_tether_lengths(integrator)
    s.reel_out_speeds  .= s.get_tether_speeds(integrator)
    s.L_C               = s.get_L_C(integrator)
    s.L_D               = s.get_L_D(integrator)
    s.D_C               = s.get_D_C(integrator)
    s.D_D               = s.get_D_D(integrator)
    calc_kite_ref_frame!(s, s.pos[s.num_E], s.pos[s.num_C], s.pos[s.num_D])

    @assert all(abs.(s.steering_pos) .<= s.set.tip_length)
    nothing
end

function model!(s::KPS4_3L, pos_; torque_control=false)
    pos2_ = zeros(3, s.num_A)
    [pos2_[:,i] .= pos_[i] for i in 1:s.num_A]
    @parameters begin
        set_values[1:3] = s.set_values
        v_wind_gnd[1:3] = s.v_wind_gnd
        v_wind[1:3] = s.v_wind
        stiffness_factor = s.stiffness_factor
    end
    @variables begin
        pos(t)[1:3, 1:s.num_A] = pos2_ # left right middle
        vel(t)[1:3, 1:s.num_A] = zeros(3, s.num_A) # left right middle
        acc(t)[1:3, 1:s.num_A] = zeros(3, s.num_A) # left right middle
        tether_length(t)[1:3]  = s.tether_lengths
        steering_pos(t)[1:2]   = s.steering_pos
        steering_vel(t)[1:2]   = zeros(2)
        steering_acc(t)[1:2]   = zeros(2)
        tether_speed(t)[1:3]   = zeros(3) # left right middle
        segment_length(t)[1:3] = zeros(3) # left right middle
        mass_tether_particle(t)[1:3]      # left right middle
        damping(t)[1:3] = s.set.damping ./ s.tether_lengths ./ s.set.segments   # left right middle
        c_spring(t)[1:3] = s.set.c_spring ./ s.tether_lengths ./ s.set.segments # left right middle
        P_c(t)[1:3] = 0.5 .* (s.pos[s.num_C] + s.pos[s.num_D])
        e_x(t)[1:3]
        e_y(t)[1:3]
        e_z(t)[1:3]
        force(t)[1:3, 1:s.num_A] = zeros(3, s.num_A) # left right middle
        rho_kite(t) = 0.0
    end
    # Collect the arrays into variables
    pos = collect(pos)
    vel = collect(vel)
    acc = collect(acc)

    eqs1 = []
    mass_per_meter = s.set.rho_tether * π * (s.set.d_tether/2000.0)^2

    [eqs1 = vcat(eqs1, D.(pos[:, i]) .~ 0.0) for i in 1:3]
    [eqs1 = vcat(eqs1, D.(pos[:, i]) .~ vel[:,i]) for i in 4:s.num_E-3]
    eqs1 = [eqs1; D.(steering_pos)   .~ steering_vel]
    [eqs1 = vcat(eqs1, D.(pos[:, i]) .~ vel[:,i]) for i in s.num_E:s.num_A]
    [eqs1 = vcat(eqs1, D.(vel[:, i]) .~ 0.0) for i in 1:3]
    [eqs1 = vcat(eqs1, D.(vel[:, i]) .~ acc[:,i]) for i in 4:s.num_E-3]
    eqs1 = [eqs1; D.(steering_vel)   .~ steering_acc]
    [eqs1 = vcat(eqs1, D.(vel[:, i]) .~ acc[:,i]) for i in s.num_E:s.num_A]

    eqs1 = vcat(eqs1, D.(tether_length) .~ tether_speed)
    if torque_control
        eqs1 = vcat(eqs1, D.(tether_speed) .~ [calc_acc_torque(tether_speed[i], norm(force[:, (i-1) % 3 + 1]),
                                                               set_values[i]) for i in 1:3])
    else
        eqs1 = vcat(eqs1, D.(tether_speed) .~ [calc_acc_speed(tether_speed[i], norm(force[:,(i-1) % 3 + 1]), 
                                                               set_values[i]) for i in 1:3])
    end

    # Compute the masses and forces
    force_eqs = SizedArray{Tuple{3, s.num_A}, Symbolics.Equation}(undef)
    force_eqs[:, :] .= (force[:, :] .~ 0)

    eqs2 = [
        pos[:, s.num_E-2] ~ pos[:, s.num_C] + e_z * steering_pos[1]
        pos[:, s.num_E-1] ~ pos[:, s.num_D] + e_z * steering_pos[2]
        vel[:, s.num_E-2] ~ vel[:, s.num_C] + e_z * steering_vel[1]
        vel[:, s.num_E-1] ~ vel[:, s.num_D] + e_z * steering_vel[2]
        acc[:, s.num_E-2] ~ acc[:, s.num_C] + e_z * steering_acc[1]
        acc[:, s.num_E-1] ~ acc[:, s.num_D] + e_z * steering_acc[2]
        segment_length       ~ tether_length  ./ s.set.segments
        mass_tether_particle ~ mass_per_meter .* segment_length
        damping              ~ s.set.damping  ./ segment_length
        c_spring             ~ s.set.c_spring ./ segment_length
        P_c ~ 0.5 * (pos[:, s.num_C] + pos[:, s.num_D])
        e_y ~ (pos[:, s.num_C] - pos[:, s.num_D]) / norm(pos[:, s.num_C] - pos[:, s.num_D])
        e_z ~ (pos[:, s.num_E] - P_c) / norm(pos[:, s.num_E] - P_c)
        e_x ~ cross(e_y, e_z)
        rho_kite ~ calc_rho(s.am, pos[3,s.num_A])
    ]

    eqs2, force_eqs = calc_aero_forces_mtk!(s, eqs2, force_eqs, force, pos, vel, t, e_x, e_y, e_z, rho_kite, v_wind)
    eqs2, force_eqs = inner_loop_mtk!(s, eqs2, force_eqs, t, force, pos, vel, segment_length, c_spring, damping, 
                                      v_wind_gnd, stiffness_factor)
    
    for i in 1:3
        eqs2 = vcat(eqs2, vcat(force_eqs[:, i]))
        eqs2 = vcat(eqs2, acc[:, i] .~ 0)
    end
    for i in 4:s.num_E-3
        eqs2 = vcat(eqs2, vcat(force_eqs[:,i]))
        eqs2 = vcat(eqs2, acc[:, i] .~ [0.0; 0.0; -G_EARTH] .+ (force[:,i] ./ mass_tether_particle[(i-1)%3+1]))
    end
    for i in s.num_E-2:s.num_E-1
        flap_resistance   = [50.0 * ((vel[:,i]-vel[:, s.num_C]) ⋅ e_z) * e_z[j] for j in 1:3]
        [force_eqs[j,i]   = force[j,i] ~ force_eqs[j,i].rhs + [0.0; 0.0; -G_EARTH][j] + flap_resistance[j] for j in 1:3]
        tether_rhs        = [force_eqs[j, i].rhs for j in 1:3]
        kite_rhs          = [force_eqs[j, i+3].rhs for j in 1:3]
        f_xy              = (tether_rhs ⋅ e_z) .* e_z
        force_eqs[:,i]   .= force[:, i] .~ tether_rhs .- f_xy
        force_eqs[:,i+3] .= force[:, i+3] .~ kite_rhs .+ f_xy
        eqs2              = vcat(eqs2, vcat(force_eqs[:, i]))
        eqs2              = vcat(eqs2, steering_acc[i-s.num_E+3] ~ (force[:,i] ./ mass_tether_particle[(i-1) % 3 + 1]) ⋅ 
                                                                    e_z - (acc[:, i+3] ⋅ e_z))
    end
    for i in s.num_E:s.num_A
        eqs2 = vcat(eqs2, vcat(force_eqs[:,i]))
        eqs2 = vcat(eqs2, acc[:, i] .~ [0.0; 0.0; -G_EARTH] .+ (force[:, i] ./ s.masses[i]))
    end

    eqs = vcat(eqs1, eqs2)

    println("making mtk model")
    @time @named sys = ODESystem(Symbolics.scalarize.(reduce(vcat, Symbolics.scalarize.(eqs))), t)
    println("making simple sys")
    @time simple_sys = structural_simplify(sys)
    return simple_sys, sys
end

# function steady_state_model!(s::KPS4_3L, pos_)
#     pos_xz_ = zeros(2, div(s.num_A, 3) * 2)
#     pos_y_ = zeros(s.num_A)
#     # [pos_xz_[:,i] .= [pos_[i][1], pos_[i][3]] for i in 1:s.num_A if i%3 == 1 || i%3 == 0]
#     j = 1
#     for i in 1:s.num_A
#         if i%3 == 1 || i%3 == 0
#             pos_xz_[:,j] .= [pos_[i][1], pos_[i][3]]
#             println(pos_xz_[:,j])
#             j += 1
#         end
#     end
#     [pos_y_[i] = pos_[i][2] for i in 1:s.num_A]
#     @parameters begin
#         set_values[1:3] = s.set_values
#         v_wind_gnd[1:3] = s.v_wind_gnd
#         v_wind[1:3] = s.v_wind
#         stiffness_factor = s.stiffness_factor
#         pos_y[1:s.num_A] = pos_y_
#     end
#     @variables begin
#         pos_xz(t)[1:2, 1:div(s.num_A, 3) * 2] = pos_xz_ # [x_left, z_left], [x_middle, z_middle]
#         vel(t)[1:3, 1:s.num_A] = zeros(3, s.num_A) # left right middle
#         acc(t)[1:3, 1:s.num_A] = zeros(3, s.num_A) # left right middle
#         tether_length(t)[1:3]  = s.tether_lengths
#         steering_pos(t)   = s.steering_pos[1]
#         steering_vel(t)   = 0.0
#         steering_acc(t)   = 0.0
#         tether_speed(t)[1:3]   = zeros(3) # left right middle
#         segment_length(t)[1:3] = zeros(3) # left right middle
#         mass_tether_particle(t)[1:3]      # left right middle
#         damping(t)[1:3] = s.set.damping ./ s.tether_lengths ./ s.set.segments   # left right middle
#         c_spring(t)[1:3] = s.set.c_spring ./ s.tether_lengths ./ s.set.segments # left right middle
#         P_c(t)[1:3] = 0.5 .* (s.pos[s.num_C] + s.pos[s.num_D])
#         e_x(t)[1:3]
#         e_y(t)[1:3]
#         e_z(t)[1:3]
#         force(t)[1:3, 1:s.num_A] = zeros(3, s.num_A) # left right middle
#         rho_kite(t) = 0.0
#     end
#     # Collect the arrays into variables
#     pos_xz = collect(pos_xz)
#     vel = collect(vel)
#     acc = collect(acc)

#     pos = Array{Union{Float64, Symbolics.Num}}(undef, 3, s.num_A)

#     # [pos[:,i] .= [pos_xz[1,i], pos_y[i], pos_xz[2,i]] for i in 1:s.num_A]
#     j = 1
#     for i in 1:s.num_A
#         if i%3 == 1
#             pos[:,i] .= [pos_xz[1,j], pos_y[i], pos_xz[2,j]] # left tether
#         elseif i%3 == 2
#             pos[:,i] .= [pos_xz[1,j], pos_y[i], pos_xz[2,j]] # right tether == -left tether
#             j += 1
#             println(pos[:,i])
#         elseif i%3 == 0
#             pos[:,i] .= [pos_xz[1,j], pos_y[i], pos_xz[2,j]] # middle tether
#             j += 1
#             println(pos[:,i])
#         end
#     end

#     eqs1 = []
#     mass_per_meter = s.set.rho_tether * π * (s.set.d_tether/2000.0)^2

#     [eqs1 = vcat(eqs1, D.(pos[1, i]) ~ 0.0) for i in 1:3 if i%3 == 1 || i%3 == 0]
#     [eqs1 = vcat(eqs1, D.(pos[3, i]) ~ 0.0) for i in 1:3 if i%3 == 1 || i%3 == 0]
#     [eqs1 = vcat(eqs1, D.(pos[1, i]) ~ vel[1,i]) for i in 4:s.num_E-3 if i%3 == 1 || i%3 == 0]
#     [eqs1 = vcat(eqs1, D.(pos[3, i]) ~ vel[3,i]) for i in 4:s.num_E-3 if i%3 == 1 || i%3 == 0]
#     eqs1 = [eqs1; D.(steering_pos)   .~ steering_vel]
#     [eqs1 = vcat(eqs1, D.(pos[1, i]) ~ vel[1,i]) for i in s.num_E:s.num_A if i%3 == 1 || i%3 == 0]
#     [eqs1 = vcat(eqs1, D.(pos[3, i]) ~ vel[3,i]) for i in s.num_E:s.num_A if i%3 == 1 || i%3 == 0]
#     [eqs1 = vcat(eqs1, D.(vel[:, i]) .~ 0.0) for i in 1:3]
#     [eqs1 = vcat(eqs1, D.(vel[:, i]) .~ acc[:,i]) for i in 4:s.num_E-3]
#     eqs1 = [eqs1; D.(steering_vel)   .~ steering_acc]
#     [eqs1 = vcat(eqs1, D.(vel[:, i]) .~ acc[:,i]) for i in s.num_E:s.num_A]

#     eqs1 = vcat(eqs1, D.(tether_length) .~ tether_speed)
#     eqs1 = vcat(eqs1, D.(tether_speed) .~ 0.0)

#     # Compute the masses and forces
#     force_eqs = SizedArray{Tuple{3, s.num_A}, Symbolics.Equation}(undef)
#     force_eqs[:, :] .= (force[:, :] .~ 0)

#     vel[:, s.num_E-2] ~ vel[:, s.num_C] + e_z * steering_vel
#     eqs2 = [
#         pos[1, s.num_E-2]   ~ pos[1, s.num_C] + e_z[1] * steering_pos
#         pos[3, s.num_E-2]   ~ pos[3, s.num_C] + e_z[3] * steering_pos
#         vel[:, s.num_E-2]   ~ vel[:, s.num_C] + e_z * steering_vel
#         vel[:, s.num_E-1]   ~ vel[:, s.num_D] + e_z * steering_vel
#         acc[:, s.num_E-2]   ~ acc[:, s.num_C] + e_z * steering_acc
#         acc[:, s.num_E-1]   ~ acc[:, s.num_D] + e_z * steering_acc
#         segment_length       ~ tether_length  ./ s.set.segments
#         mass_tether_particle ~ mass_per_meter .* segment_length
#         damping              ~ s.set.damping  ./ segment_length
#         c_spring             ~ s.set.c_spring ./ segment_length
#         P_c                 ~ 0.5 * (pos[:, s.num_C] + pos[:, s.num_D])
#         e_y                 ~ (pos[:, s.num_C] - pos[:, s.num_D]) / norm(pos[:, s.num_C] - pos[:, s.num_D])
#         e_z                 ~ (pos[:, s.num_E] - P_c) / norm(pos[:, s.num_E] - P_c)
#         e_x                 ~ cross(e_y, e_z)
#         rho_kite ~ calc_rho(s.am, pos[3,s.num_A])
#     ]

#     eqs2, force_eqs = calc_aero_forces_mtk!(s, eqs2, force_eqs, force, pos, vel, t, e_x, e_y, e_z, rho_kite, v_wind)
#     eqs2, force_eqs = inner_loop_mtk!(s, eqs2, force_eqs, t, force, pos, vel, segment_length, c_spring, damping, 
#                                       v_wind_gnd, stiffness_factor)
    
#     for i in 1:3
#         eqs2 = vcat(eqs2, vcat(force_eqs[:, i]))
#         eqs2 = vcat(eqs2, acc[:, i] .~ 0)
#     end
#     for i in 4:s.num_E-3
#         eqs2 = vcat(eqs2, vcat(force_eqs[:,i]))
#         eqs2 = vcat(eqs2, acc[:, i] .~ [0.0; 0.0; -G_EARTH] .+ (force[:,i] ./ mass_tether_particle[(i-1)%3+1]))
#     end
#     for i in s.num_E-2:s.num_E-1
#         flap_resistance   = [50.0 * ((vel[:,i]-vel[:, s.num_C]) ⋅ e_z) * e_z[j] for j in 1:3]
#         [force_eqs[j,i]   = force[j,i] ~ force_eqs[j,i].rhs + [0.0; 0.0; -G_EARTH][j] + flap_resistance[j] for j in 1:3]
#         tether_rhs        = [force_eqs[j, i].rhs for j in 1:3]
#         kite_rhs          = [force_eqs[j, i+3].rhs for j in 1:3]
#         f_xy              = (tether_rhs ⋅ e_z) .* e_z
#         force_eqs[:,i]   .= force[:, i] .~ tether_rhs .- f_xy
#         force_eqs[:,i+3] .= force[:, i+3] .~ kite_rhs .+ f_xy
#         eqs2              = vcat(eqs2, vcat(force_eqs[:, i]))
#     end
#     eqs2              = vcat(eqs2, steering_acc ~ (force[:,s.num_E-2] ./ mass_tether_particle[1]) ⋅ 
#                                                                 e_z - (acc[:, s.num_C] ⋅ e_z))
#     for i in s.num_E:s.num_A
#         eqs2 = vcat(eqs2, vcat(force_eqs[:,i]))
#         eqs2 = vcat(eqs2, acc[:, i] .~ [0.0; 0.0; -G_EARTH] .+ (force[:, i] ./ s.masses[i]))
#     end

#     eqs = vcat(eqs1, eqs2)

#     println("making mtk model")
#     @time @named sys = ODESystem(Symbolics.scalarize.(reduce(vcat, Symbolics.scalarize.(eqs))), t)
#     println("making simple sys")
#     @time simple_sys = structural_simplify(sys)
#     return simple_sys, sys
# end