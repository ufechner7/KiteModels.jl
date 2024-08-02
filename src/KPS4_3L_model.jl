


"""
    calc_aero_forces!(s::KPS4_3L, pos, vel)

Calculates the aerodynamic forces acting on the kite particles.

Parameters:
- pos:              vector of the particle positions
- vel:              vector of the particle velocities
- rho:              air density [kg/m^3]
- rel_depower:      value between  0.0 and  1.0
- alpha_depower:    depower angle [degrees]
- rel_steering:     value between -1.0 and +1.0

Updates the vector s.forces of the first parameter.
"""
function calc_aero_forces_model!(s::KPS4_3L, eqs2, force_eqs, force, pos, vel, t, e_x, e_y, e_z)
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
    E_c = collect(E_c)
    v_cx = collect(v_cx)
    v_dx = collect(v_dx)
    v_dy = collect(v_dy)
    v_dz = collect(v_dz)
    v_cy = collect(v_cy)
    v_cz = collect(v_cz)

    eqs2 = [
        eqs2
        δ_left ~ (pos[:,s.num_E-2].-pos[:,s.num_C]) ⋅ e_z
        δ_right ~ (pos[:,s.num_E-1].-pos[:,s.num_D]) ⋅ e_z
        E_c .~ pos[:,s.num_E] .+ e_z .* (-s.set.bridle_center_distance + s.set.radius) # in the aero calculations, E_c is the center of the circle shape on which the kite lies
        v_cx .~ dot(vel[:,s.num_C], e_x).*e_x
        v_dx .~ dot(vel[:,s.num_D], e_x).*e_x
        v_dy .~ dot(vel[:,s.num_D], e_y).*e_y
        v_dz .~ dot(vel[:,s.num_D], e_z).*e_z
        v_cy .~ dot(vel[:,s.num_C], e_y).*e_y
        v_cz .~ dot(vel[:,s.num_C], e_z).*e_z
        y_lc ~ norm(pos[:,s.num_C] .- 0.5 .* (pos[:,s.num_C].+pos[:,s.num_D]))
        y_ld ~ -norm(pos[:,s.num_D] .- 0.5 .* (pos[:,s.num_C].+pos[:,s.num_D]))
    ]

    # integrating loop variables
    @variables begin
        F(t)[1:3, 1:n*2]
        e_r(t)[1:3, 1:n*2]
        y_l(t)[1:n*2]
        v_kite(t)[1:3, 1:n*2]
        v_a(t)[1:3, 1:n*2]
        e_drift(t)[1:3, 1:n*2]
        v_a_xr(t)[1:3, 1:n*2]
        aoa(t)[1:n*2]
        dL_dα(t)[1:3, 1:n*2]
        dD_dα(t)[1:3, 1:n*2]
        L_C(t)[1:3]
        L_D(t)[1:3]
        D_C(t)[1:3]
        D_D(t)[1:3]
        F_steering_c(t)[1:3]
        F_steering_d(t)[1:3]
        d(t)[1:n*2]
    end
    F = collect(F)
    e_r = collect(e_r)
    y_l = collect(y_l)
    v_kite = collect(v_kite)
    v_a = collect(v_a)
    e_drift = collect(e_drift)
    v_a_xr = collect(v_a_xr)
    aoa = collect(aoa)
    dL_dα = collect(dL_dα)
    dD_dα = collect(dD_dα)
    l_c_eq = collect(L_C .~ 0)
    l_d_eq = collect(L_D .~ 0)
    d_c_eq = collect(D_C .~ 0)
    d_d_eq = collect(D_D .~ 0)
    F_steering_c = collect(F_steering_c)
    F_steering_d = collect(F_steering_d)
    kite_length = zeros(MVector{n*2, SimFloat})
    α = zero(SimFloat)
    α_0 = zero(SimFloat)
    α_middle = zero(SimFloat)
    dα = zero(SimFloat)
    # Calculate the lift and drag
    α_0 ~ pi/2 - s.set.width/2/s.set.radius
    α_middle ~ pi/2
    dα ~ (α_middle - α_0) / n
    # println("calculating aero forces...")
    for i in 1:n*2
        if i <= n
            α = α_0 + -dα/2 + i*dα
        else
            α = pi - (α_0 + -dα/2 + (i-n)*dα)
        end

        eqs2 = [
            eqs2
            F[:,i] ~ E_c .+ e_y.*cos(α).*s.set.radius .- e_z.*sin(α).*s.set.radius
            e_r[:,i] ~ (E_c .- F[:,i])./norm(E_c .- F[:,i])
            y_l[i] ~ cos(α) * s.set.radius
            α < π/2 ?
                v_kite[:,i] ~ ((v_cx .- v_dx)./(y_lc .- y_ld).*(y_l[i] .- y_ld) .+ v_dx) .+ v_cy .+ v_cz :
                v_kite[:,i] ~ ((v_cx .- v_dx)./(y_lc .- y_ld).*(y_l[i] .- y_ld) .+ v_dx) .+ v_dy .+ v_dz
            v_a[:,i] ~ s.v_wind .- v_kite[:,i]
            e_drift[:,i] ~ (e_r[:,i] × e_x)
            v_a_xr[:,i] ~ v_a[:,i] .- (v_a[:,i] ⋅ e_drift[:,i]) .* e_drift[:,i]
        ]
        if α < π/2
            kite_length[i] = (s.set.tip_length + (s.set.middle_length-s.set.tip_length)*α*s.set.radius/(0.5*s.set.width))
        else
            kite_length[i] = (s.set.tip_length + (s.set.middle_length-s.set.tip_length)*(π-α)*s.set.radius/(0.5*s.set.width))
        end
        eqs2 = [
            eqs2
            α < s.α_l ?
                d[i] ~ δ_left :
            α > s.α_r ?
                d[i] ~ δ_right :
                d[i] ~ (δ_right - δ_left) / (α_r - α_l) * (α - α_l) + (δ_left)
            aoa[i] ~ π - acos2((v_a_xr[:,i]./norm(v_a_xr[:,i])) ⋅ e_x) + asin(clamp(d[i]/kite_length[i], -1.0, 1.0))
            dL_dα[:,i] .~ 0.5*s.rho*(norm(v_a_xr[:,i]))^2*s.set.radius*kite_length[i]*rad_cl_model(aoa[i]) .* ((v_a_xr[:,i] × e_drift[:,i]) ./ norm(v_a_xr[:,i] × e_drift[:,i]))
            dD_dα[:,i] .~ 0.5*s.rho*norm(v_a_xr[:,i])*s.set.radius*kite_length[i]*rad_cd_model(aoa[i]) .* v_a_xr[:,i] # the sideways drag cannot be calculated with the C_d formula
        ]
        if i <= n
            [l_c_eq[j] = (L_C[j] ~ l_c_eq[j].rhs + dL_dα[j,i] * dα) for j in 1:3]
            [d_c_eq[j] = (D_C[j] ~ d_c_eq[j].rhs + dD_dα[j,i] * dα) for j in 1:3]
        else 
            [l_d_eq[j] = (L_D[j] ~ l_d_eq[j].rhs + dL_dα[j,i] * dα) for j in 1:3]
            [d_d_eq[j] = (D_D[j] ~ d_d_eq[j].rhs + dD_dα[j,i] * dα) for j in 1:3]
        end
    end
    
    eqs2 = [
        eqs2
        l_c_eq
        d_c_eq
        l_d_eq
        d_d_eq
        F_steering_c ~ ((0.1 * (L_C ⋅ -e_z)) .* -e_z)
        F_steering_d ~ ((0.1 * (L_D ⋅ -e_z)) .* -e_z)
        ]
    
    force_eqs[:,s.num_C] .= (force[:,s.num_C] .~ (L_C + D_C) .- F_steering_c) 
    force_eqs[:,s.num_D] .= (force[:,s.num_D] .~ (L_D .+ D_D) .- F_steering_d)
    force_eqs[:,s.num_E-2] .= (force[:,s.num_E-2] .~ F_steering_c)
    force_eqs[:,s.num_E-1] .= (force[:,s.num_E-1] .~ F_steering_d)
    return eqs2, force_eqs
end

""" 
    calc_particle_forces!(s::KPS4_3L, pos1, pos2, vel1, vel2, spring, d_tether, rho, i)

Calculate the drag force and spring force of the tether segment, defined by the parameters pos1, pos2, vel1 and vel2
and distribute it equally on the two particles, that are attached to the segment.
The result is stored in the array s.forces. 
"""
@inline function calc_particle_forces_model!(s::KPS4_3L, eqs2, force_eqs, force, pos1, pos2, vel1, vel2, length, c_spring, damping, d_tether, rho, i,
    l_0, k, c, segment, rel_vel, av_vel, norm1, unit_vector, k1, k2, c1, spring_vel,
            spring_force, v_apparent, area, v_app_perp, half_drag_force)
    eqs2 = [
        eqs2
        i <= s.set.segments*3 ? l_0 ~ length[i%3+1] : l_0 ~ s.springs[i].length # Unstressed length
        i <= s.set.segments*3 ? k ~ c_spring[i%3+1] * s.stiffness_factor : k ~ s.springs[i].c_spring * s.stiffness_factor # Spring constant
        i <= s.set.segments*3 ? c ~ damping[i%3+1] : c ~ s.springs[i].damping                    # Damping coefficient    
        segment .~ pos1 - pos2
        rel_vel .~ vel1 - vel2
        av_vel .~ 0.5 * (vel1 + vel2)
        norm1 ~ norm(segment)
        unit_vector .~ segment / norm1
        k1 ~ 0.25 * k # compression stiffness kite segments
        k2 ~ 0.1 * k  # compression stiffness tether segments
        c1 ~ 6.0 * c  # damping kite segments
        spring_vel .~ rel_vel ⋅ unit_vector
    ]

    if i > s.num_E  # kite springs
        for j in 1:3
            eqs2 = [
                eqs2
                spring_force[j] ~ ifelse(
                    (norm1 - l_0) > 0.0,
                    (k*(l_0 - norm1) - c1 * spring_vel) * unit_vector[j],
                    (k1*(l_0 - norm1) - c * spring_vel) * unit_vector[j]
                )
            ]
        end
    else
        for j in 1:3
            eqs2 = [
                eqs2
                spring_force[j] ~ ifelse(
                    (norm1 - l_0) > 0.0,
                    (k*(l_0 - norm1) - c * spring_vel) * unit_vector[j],
                    (k2*(l_0 - norm1) - c * spring_vel) * unit_vector[j]
                    )
            ]
        end
    end
    eqs2 = [
        eqs2
        v_apparent .~ s.v_wind_tether .- av_vel
        area ~ norm1 * d_tether
        v_app_perp ~ s.v_apparent .- s.v_apparent ⋅ unit_vector * unit_vector
        half_drag_force .~ (0.25 * rho * s.set.cd_tether * norm(v_app_perp) * area) .* v_app_perp
    ]

    for j in 1:3
        force_eqs[j, s.springs[i].p1] = (force[j,s.springs[i].p1] ~ force_eqs[j, s.springs[i].p1].rhs + half_drag_force[j] + spring_force[j])
        force_eqs[j, s.springs[i].p2] = (force[j,s.springs[i].p2] ~ force_eqs[j, s.springs[i].p2].rhs + half_drag_force[j] - spring_force[j])
    end
    
    # if i <= 3 @inbounds s.last_forces[i%3+1] .~ s.forces[spring.p1] end
    return eqs2, force_eqs
end

"""calc_aero_forces
    inner_loop_model!(s::KPS4_3L, pos, vel, v_wind_gnd, segments, d_tether)

Calculate the forces, acting on all particles.

Output:length
- s.forces
- s.v_wind_tether
"""
@inline function inner_loop_model!(s::KPS4_3L, eqs2, force_eqs, t, force, pos, vel, length, c_spring, damping, v_wind_gnd, d_tether)
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
    v_wind_tether = collect(v_wind_tether)
    segment = collect(segment)
    rel_vel = collect(rel_vel)
    av_vel = collect(av_vel)
    unit_vector = collect(unit_vector)
    spring_force = collect(spring_force)
    v_apparent = collect(v_apparent)
    v_app_perp = collect(v_app_perp)
    half_drag_force = collect(half_drag_force)
    
    for i in eachindex(s.springs)
        p1 = s.springs[i].p1  # First point nr.
        p2 = s.springs[i].p2
        eqs2 = [
            eqs2
            height[i] ~ 0.5 * (pos[:,p1][3] + pos[:,p2][3])
            rho[i] ~ calc_rho(s.am, height[i])
            v_wind_tether[:,i] .~ calc_wind_factor(s.am, height[i]) .* v_wind_gnd
        ]
        # TODO: @assert height > 0
        eqs2, force_eqs = calc_particle_forces_model!(s, eqs2, force_eqs, force, pos[:,p1], pos[:,p2], vel[:,p1], vel[:,p2], length, c_spring, damping, d_tether, rho[i], i,
            l_0[i], k[i], c[i], segment[:,i], rel_vel[:,i], av_vel[:,i], norm1[i], 
            unit_vector[:,i], k1[i], k2[i], c1[i], spring_vel[i],
            spring_force[:,i], v_apparent[:,i], area[i], v_app_perp[:,i], half_drag_force[:,i])
    end

    return eqs2, force_eqs
end

function model!(s::KPS4_3L, pos_, vel_)
    pos2_ = zeros(3, s.num_A)
    vel2_ = zeros(3, s.num_A)
    [pos2_[:,i] .= pos_[i] for i in 1:s.num_A]
    [vel2_[:,i] .= vel_[i] for i in 1:s.num_A]
    println(size(pos2_))
    @independent_variables t
    @variables begin
        pos(t)[1:3, 1:s.num_A] = pos2_
        vel(t)[1:3, 1:s.num_A] = vel2_
        acc(t)[1:3, 1:s.num_A] = zeros(3, s.num_A)
        lengths(t)[1:3] = ones(3)
        steering_pos(t)[1:2] = zeros(2)
        steering_vel(t)[1:2] = ones(2)
        steering_acc(t)[1:2] = ones(2)
        reel_out_speed(t)[1:3] = ones(3)
        segment_lengths(t)[1:3] = ones(3)
        mass_tether_particle(t)[1:3] = ones(3)
        damping(t)[1:3] = ones(3)
        c_spring(t)[1:3] = ones(3)
        P_c(t)[1:3] = ones(3)
        e_x(t)[1:3] = ones(3)
        e_y(t)[1:3] = ones(3)
        e_z(t)[1:3] = ones(3)
        force(t)[1:3, 1:s.num_A] = ones(3, s.num_A)
    end
    # Collect the arrays into variables
    pos = collect(pos)
    vel = collect(vel)
    acc = collect(acc)
    lengths = collect(lengths)
    steering_pos = collect(steering_pos)
    steering_vel = collect(steering_vel)
    steering_acc = collect(steering_acc)
    reel_out_speed = collect(reel_out_speed)
    segment_lengths = collect(segment_lengths)
    mass_tether_particle = collect(mass_tether_particle)
    damping = collect(damping)
    c_spring = collect(c_spring)
    P_c = collect(P_c)
    e_x = collect(e_x)
    e_y = collect(e_y)
    e_z = collect(e_z)
    force = collect(force)

    D = Differential(t)

    eqs1 = []
    mass_per_meter = s.set.rho_tether * π * (s.set.d_tether/2000.0)^2

    for i in 1:3
        eqs1 = vcat(
            eqs1,
            D.(pos[:,i]) .~ 0.0, # fix pos of s.num_E-2 and s.num_E-1
            D.(vel[:,i]) .~ 0.0
        )
    end
    for i in 4:s.num_E-3
        eqs1 = vcat(
            eqs1,
            D.(pos[:,i]) .~ vel[:,i],
            D.(vel[:,i]) .~ acc[:,i]
        )
    end
    eqs1 = [
        eqs1
        D.(steering_pos) .~ steering_vel
        D.(steering_vel) .~ steering_acc
        pos[:,s.num_E-2] .~ pos[:,s.num_C] .+ e_z .* steering_pos[1]
        pos[:,s.num_E-1] .~ pos[:,s.num_D] .+ e_z .* steering_pos[2]
        vel[:,s.num_E-2] .~ vel[:,s.num_C] .+ e_z .* steering_vel[1]
        vel[:,s.num_E-1] .~ vel[:,s.num_D] .+ e_z .* steering_vel[2]
        acc[:,s.num_E-2] .~ acc[:,s.num_C] .+ e_z .* steering_acc[1]
        acc[:,s.num_E-1] .~ acc[:,s.num_D] .+ e_z .* steering_acc[2]
    ]
    for i in s.num_E:s.num_A
        eqs1 = vcat(
            eqs1,
            D.(pos[:,i]) .~ vel[:,i],
            D.(vel[:,i]) .~ acc[:,i]
        )
    end
    eqs1 = [
        eqs1
        D.(lengths) .~ reel_out_speed
        D.(reel_out_speed) .~ 0 # TODO: add winch function
    ]

    # Compute the masses and forces
    eqs2 = []
    force_eqs::SizedArray{Tuple{3, s.num_A}, Symbolics.Equation} = SizedArray{Tuple{3, s.num_A}, Symbolics.Equation}(undef)
    force_eqs[:,:] .= (force[:,:] .~ 0)

    eqs2 = [
        eqs2
        segment_lengths .~ lengths ./ s.set.segments
        mass_tether_particle .~ mass_per_meter .* segment_lengths
        damping .~ s.set.damping ./ segment_lengths
        c_spring .~ s.set.c_spring ./ segment_lengths
        P_c ~ 0.5 .* (pos[:,s.num_C]+pos[:,s.num_D])
        # e_y .~ norm(P_c)
        e_y .~ (pos[:,s.num_C] .- pos[:,s.num_D]) ./ norm(pos[:,s.num_C] .- pos[:,s.num_D])
        e_z .~ (pos[:,s.num_E] .- P_c) ./ norm(pos[:,s.num_E] .- P_c)
        e_x .~ cross(e_y, e_z)
    ]

    eqs2, force_eqs = calc_aero_forces_model!(s, eqs2, force_eqs, force, pos, vel, t, e_x, e_y, e_z)
    eqs2, force_eqs = inner_loop_model!(s, eqs2, force_eqs, t, force, pos, vel, s.v_wind_gnd, segment_lengths, c_spring, damping, s.set.d_tether/1000.0)
    
    for i in 1:3
        eqs2 = vcat(eqs2, acc[:,i] .~ 0)
    end
    for i in s.num_E-2:s.num_E-1
        # println(force[1,i].rhs)
        [force_eqs[j,i] = force[j,i] ~ force_eqs[j,i].rhs + [0, 0, -G_EARTH][j] + 500.0 * ((vel[:,i]-vel[:,s.num_C]) ⋅ e_z) * e_z[j] for j in 1:3] # TODO: more damping
        [force_eqs[j,i] = force[j,i] ~ force_eqs[j,i].rhs - s.forces[i][j] - (s.forces[i] ⋅ e_z) * e_z[j] for j in 1:3]
        [force_eqs[j,i+3] = force[j,i+3] ~ force_eqs[j,i].rhs + s.forces[i][j] - (s.forces[i] ⋅ e_z) * e_z[j] for j in 1:3]
        eqs2 = vcat(eqs2, vcat(force_eqs[:,i]))
        # println(size((force[:,i] ./ mass_tether_particle[i%3+1])) ⋅ e_z)
        eqs2 = vcat(eqs2, steering_acc[i-s.num_E+3] ~ (force[:,i] ./ mass_tether_particle[i%3+1]) ⋅ e_z - (acc[:,i+3] ⋅ s.e_z))
    end
    for i in 4:s.num_E-3
        eqs2 = vcat(eqs2, vcat(force_eqs[:,i]))
        eqs2 = vcat(eqs2, acc[:,i] .~ (force[:,i] ./ mass_tether_particle[i%3+1]))
    end
    for i in s.num_E:s.num_A
        eqs2 = vcat(eqs2, vcat(force_eqs[:,i]))
        eqs2 = vcat(eqs2, acc[:,i] .~ (force[:,i] ./ s.masses[i]))
    end

    eqs = vcat(eqs1, eqs2)

    # @assert false
    println("making model")
    @time @named sys = ODESystem(eqs, t)
    println("making simple sys")
    @time simple_sys = structural_simplify(sys)
    return simple_sys, sys
end
