# ==================== mtk model functions ================================================
# Implementation of the three-line model using ModellingToolkit.jl

function calc_speed_acc(winch::AsyncMachine, tether_vel, norm_, set_speed)
    calc_acceleration(winch, tether_vel, norm_; set_speed, set_torque=nothing, use_brake=false) # TODO: add brake setting
end
function calc_torque_acc(winch::TorqueControlledMachine, tether_vel, norm_, set_torque)
    calc_acceleration(winch, tether_vel, norm_; set_speed=nothing, set_torque, use_brake=false)
end
@register_symbolic calc_speed_acc(winch::AsyncMachine, tether_vel, norm_, set_speed)
@register_symbolic calc_torque_acc(winch::TorqueControlledMachine, tether_vel, norm_, set_torque)

function sym_interp(interp::Function, aoa, trailing_edge_angle)
    return interp(rad2deg(aoa), rad2deg(trailing_edge_angle-aoa)) # TODO: register callable struct https://docs.sciml.ai/Symbolics/dev/manual/functions/#Symbolics.@register_array_symbolic
end
@register_symbolic sym_interp(interp::Function, aoa, trailing_edge_angle)

function normalize(vec)
    return vec / norm(vec)
end
@register_symbolic normalize(vec)

function mean(vec::AbstractVector)
    return sum(vec) / length(vec)
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

function force_eqs!(s, system, eqs, defaults, guesses; 
        tether_kite_force, tether_kite_torque_b, R_b_w, kite_pos, q_inf, wind_vec_gnd)

    points, groups, segments, pulleys, tethers, winches = 
        system.points, system.groups, system.segments, system.pulleys, system.tethers, system.winches
    
    # ==================== POINTS ==================== #
    @variables begin
        pos(t)[1:3, eachindex(points)]
        vel(t)[1:3, eachindex(points)]
        acc(t)[1:3, eachindex(points)]
        point_force(t)[1:3, eachindex(points)]

        spring_force_vec(t)[1:3, eachindex(segments)]
        drag_force(t)[1:3, eachindex(segments)]

        twist_angle(t)[eachindex(groups)] # main body angle left / right
    end
    for point in points
        F::Vector{Num} = zeros(Num, 3)
        mass = 0.0
        for segment in segments
            if point.idx in segment.points
                mass_per_meter = s.set.rho_tether * π * (segment.diameter/2000)^2
                inverted = segment.points[2] == point.idx
                if inverted
                    F .-= spring_force_vec[:, segment.idx]
                else
                    F .+= spring_force_vec[:, segment.idx]
                end
                mass += mass_per_meter * segment.l0 / 2
                F .+= 0.5drag_force[:, segment.idx]
            end
        end
        eqs = [
            eqs
            point_force[:, point.idx]  ~ F
        ]

        winch_point = false
        for tether in tethers
            if point.idx == tether.winch_point
                winch_point = true
                break
            end
        end

        if point.type === KITE
            tether_kite_force .+= F
            tether_kite_torque_b .+= point.pos_b × (R_b_w' * F)
            chord_b = point.pos_b - point.fixed_pos

            found = 0
            group_idx = 0
            for group in groups
                if point.idx in group.points
                    group_idx = group.idx
                    found += 1
                end
            end
            !(found == 1) && throw(ArgumentError("Kite point number $(point.idx) is part of $found groups, 
                and should be part of exactly 1 groups."))

            pos_b = point.fixed_pos + rotate_v_around_k(chord_b, point.chord, twist_angle[group_idx])
            eqs = [
                eqs
                pos[:, point.idx]    ~ kite_pos + R_b_w * pos_b
                vel[:, point.idx]    ~ zeros(3)
                acc[:, point.idx]    ~ zeros(3)
            ]
        elseif point.type === WINCH
            eqs = [
                eqs
                pos[:, point.idx]    ~ point.pos_w
                vel[:, point.idx]    ~ zeros(3)
                acc[:, point.idx]    ~ zeros(3)
            ]
        elseif point.type === DYNAMIC
            eqs = [
                eqs
                D(pos[:, point.idx]) ~ vel[:, point.idx]
                D(vel[:, point.idx]) ~ acc[:, point.idx]
                acc[:, point.idx]    ~ point_force[:, point.idx] / mass + [0, 0, -G_EARTH]
            ]
            defaults = [
                defaults
                [pos[j, point.idx] => point.pos_w[j] for j in 1:3]
                [vel[j, point.idx] => 0 for j in 1:3]
            ]
        elseif point.type === STATIC
            eqs = [
                eqs
                vel[:, point.idx]    ~ zeros(3)
                acc[:, point.idx]    ~ zeros(3)
                acc[:, point.idx]    ~ point_force[:, point.idx] / mass + [0, 0, -G_EARTH]
            ]
            # guesses = [
            #     guesses
            #     [vel[j, point.idx] => 0 for j in 1:3]
            # ]
            # defaults = [
            #     defaults
            #     [pos[j, point.idx] => point.pos_w[j] for j in 1:3]
            # ]
        else
            throw(ArgumentError("Unknown point type: $(typeof(point))"))
        end
    end

    # ==================== GROUPS ==================== #
    @parameters torque_coeff_dist[eachindex(s.aero.panels)] = s.vsm_solver.sol.moment_coefficient_distribution
    @variables begin
        trailing_edge_angle(t)[eachindex(groups)] # angle left / right
        trailing_edge_ω(t)[eachindex(groups)] # angular rate
        trailing_edge_α(t)[eachindex(groups)] # angular acc
        twist_ω(t)[eachindex(groups)] # angular rate
        twist_α(t)[eachindex(groups)] # angular acc
    end
    torque_dist = torque_coeff_dist * q_inf * s.aero.projected_area
    used_panels = 0
    for group in groups
        panel_indices = Int16[]
        for (i, panel) in enumerate(s.aero.panels)
            panel_y = mean([panel.LE_point_1[2], panel.LE_point_2[2]])
            @assert group.y_lim[1] < group.y_lim[2]
            if group.y_lim[1] <= panel_y <= group.y_lim[2]
                push!(panel_indices, i)
                used_panels += 1
            end
        end
        group_torque = sum([torque_dist[i] for i in panel_indices])
        chord = points[group.points[1]].chord
        inertia = 1/3 * (s.set.mass/length(groups)) * (norm(chord))^2 # plate inertia around leading edge
        @assert !(inertia ≈ 0.0)
        eqs = [
            eqs
            D(twist_angle[group.idx]) ~ twist_ω[group.idx]
            D(twist_ω[group.idx]) ~ twist_α[group.idx]
            twist_α[group.idx] ~ group_torque / inertia
        ]
        defaults = [
            defaults
            twist_angle[group.idx] => 0
            twist_ω[group.idx] => 0
        ]
    end
    !(used_panels == length(torque_dist)) && 
        throw(ArgumentError("$used_panels out of $(length(torque_dist)) panels are used in the torque distribution"))

    # ==================== SEGMENTS ==================== #
    @variables begin
        segment_vec(t)[1:3, eachindex(segments)]
        unit_vector(t)[1:3, eachindex(segments)]
        len(t)[eachindex(segments)]
        l0(t)[eachindex(segments)]
        rel_vel(t)[1:3, eachindex(segments)]
        spring_vel(t)[eachindex(segments)]
        spring_force(t)[eachindex(segments)]
        stiffness(t)[eachindex(segments)]
        damping(t)[eachindex(segments)]

        height(t)[eachindex(segments)]
        segment_vel(t)[1:3, eachindex(segments)]
        segment_rho(t)[eachindex(segments)]
        wind_vel(t)[1:3, eachindex(segments)]
        va(t)[1:3, eachindex(segments)]
        area(t)[eachindex(segments)]
        app_perp_vel(t)[1:3, eachindex(segments)]
        drag_force(t)[1:3, eachindex(segments)]

        pulley_l0(t)[eachindex(pulleys)]

        tether_length(t)[eachindex(tethers)]
    end
    for segment in segments
        p1, p2 = segment.points[1], segment.points[2]

        guesses = [
            guesses
            [segment_vec[j, segment.idx] => points[p2].pos_w[j] - points[p1].pos_w[j] for j in 1:3]
        ]

        if segment.type === BRIDLE
            in_pulley = false
            for pulley in pulleys
                if segment.idx == pulley.segments[1] # each bridle segment has to be part of no pulley or one pulley
                    eqs = [
                        eqs
                        l0[segment.idx] ~ pulley_l0[pulley.idx]
                    ]
                    in_pulley += 1
                end
                if segment.idx == pulley.segments[2]
                    eqs = [
                        eqs
                        l0[segment.idx] ~ pulley.sum_length - pulley_l0[pulley.idx]
                    ]
                    in_pulley += 1
                end
            end
            if in_pulley == 0
                eqs = [
                    eqs
                    l0[segment.idx] ~ segment.l0
                ]
            end
            (in_pulley > 1) && throw(ArgumentError("Bridle segment number $(segment.idx) is part of
                $in_pulley pulleys, and should be part of either 0 or 1 pulleys."))
        elseif segment.type === POWER || segment.type === STEERING
            in_tether = 0
            for tether in tethers
                if segment.idx in tether.segments # each tether segment has to be part of exactly one tether
                    in_winch = 0
                    winch_idx = 0
                    for winch in winches
                        if tether.idx in winch.tethers
                            winch_idx = winch.idx
                            in_winch += 1
                        end
                    end
                    (in_winch != 1) && throw(ArgumentError("Tether number $(tether.idx) is part of
                        $(in_winch) winches, and should be part of exactly 1 winch."))

                    eqs = [
                        eqs
                        l0[segment.idx] ~ tether_length[winch_idx] / length(tether.segments)
                    ]
                    in_tether += 1
                end
            end
            (in_tether != 1) && throw(ArgumentError("Segment number $(segment.idx) is part of 
                $in_tether tethers, and should be part of exactly 1 tether."))
        else
            throw(ArgumentError("Unknown segment type: $(segment.type)"))
        end

        stiffness_m = s.set.e_tether * (segment.diameter/2000)^2 * pi
        (segment.type === BRIDLE) && (compression_frac = 1.0)
        (segment.type === POWER) && (compression_frac = 0.1)
        (segment.type === STEERING) && (compression_frac = 0.1)
        
        damping_m = (s.set.damping / s.set.c_spring) * stiffness_m
        (segment.type === BRIDLE) && (damping_m = 10damping_m)

        eqs = [
            eqs
            # spring force equations
            segment_vec[:, segment.idx]  ~ pos[:, p2] - pos[:, p1]
            len[segment.idx]             ~ norm(segment_vec[:, segment.idx])
            unit_vector[:, segment.idx]  ~ segment_vec[:, segment.idx]/len[segment.idx]
            rel_vel[:, segment.idx]      ~ vel[:, p1] - vel[:, p2]
            spring_vel[segment.idx]      ~ rel_vel[:, segment.idx] ⋅ unit_vector[:, segment.idx]
            stiffness[segment.idx]       ~ ifelse(len[segment.idx] > segment.l0,
                                        stiffness_m / len[segment.idx],
                                        compression_frac * stiffness_m / len[segment.idx]
            )
            damping[segment.idx]         ~ damping_m / len[segment.idx]
            spring_force[segment.idx] ~  (stiffness[segment.idx] * segment.l0 * (len[segment.idx] - l0[segment.idx]) - 
                            damping[segment.idx] * segment.l0 * spring_vel[segment.idx])
            spring_force_vec[:, segment.idx]  ~ spring_force[segment.idx] * unit_vector[:, segment.idx]

            # drag force equations
            height[segment.idx]          ~ max(0.0, 0.5(pos[:, p1][3] + pos[:, p2][3]))
            segment_vel[:, segment.idx]  ~ 0.5(vel[:, p1] + vel[:, p2])
            segment_rho[segment.idx]     ~ calc_rho(s.am, height[segment.idx])
            wind_vel[:, segment.idx]     ~ AtmosphericModels.calc_wind_factor(s.am, height[segment.idx], s.set.profile_law) * wind_vec_gnd
            va[:, segment.idx]           ~ wind_vel[:, segment.idx] - segment_vel[:, segment.idx]
            area[segment.idx]            ~ len[segment.idx] * segment.diameter
            app_perp_vel[:, segment.idx] ~ va[:, segment.idx] - 
                                        (va[:, segment.idx] ⋅ unit_vector[:, segment.idx]) * unit_vector[:, segment.idx]
            drag_force[:, segment.idx]   ~ (0.5 * segment_rho[segment.idx] * s.set.cd_tether * norm(va[:, segment.idx]) * 
                                        area[segment.idx]) * app_perp_vel[:, segment.idx]
        ]
    end

    # ==================== PULLEYS ==================== #
    @variables begin
        pulley_l0(t)[eachindex(pulleys)]
        pulley_vel(t)[eachindex(pulleys)]
        pulley_force(t)[eachindex(pulleys)]
        pulley_acc(t)[eachindex(pulleys)]
    end
    pulley_damping = 10
    for pulley in pulleys
        segment = segments[pulley.segments[1]]
        mass_per_meter = s.set.rho_tether * π * (segment.diameter/2000)^2
        mass = pulley.sum_length * mass_per_meter
        eqs = [
            eqs
            D(pulley_l0[pulley.idx]) ~ pulley_vel[pulley.idx]
            D(pulley_vel[pulley.idx]) ~ pulley_acc[pulley.idx]
            pulley_force[pulley.idx]    ~ spring_force[pulley.segments[1]] - spring_force[pulley.segments[2]]
            pulley_acc[pulley.idx]      ~ pulley_force[pulley.idx] / mass - pulley_damping * pulley_vel[pulley.idx]
        ]
        defaults = [
            defaults
            pulley_l0[pulley.idx] => segments[pulley.segments[1]].l0
            pulley_vel[pulley.idx] => 0
        ]
    end

    # ==================== WINCHES ==================== #
    @parameters set_values(t)[eachindex(winches)] = zeros(length(winches))
    @variables begin
        tether_vel(t)[eachindex(winches)]
        tether_acc(t)[eachindex(winches)]
        winch_force(t)[eachindex(winches)]
    end
    for winch in winches
        F = zero(Num)
        for tether_idx in winch.tethers
            point_idx = tethers[tether_idx].winch_point
            F += norm(point_force[:, point_idx])
        end
        eqs = [
            eqs
            D(tether_length[winch.idx]) ~ tether_vel[winch.idx]
            D(tether_vel[winch.idx]) ~ tether_acc[winch.idx]
            tether_acc[winch.idx] ~ calc_torque_acc( # TODO: torque and speed control
                winch.model, tether_vel[winch.idx], 
                winch_force[winch.idx], 
                set_values[winch.idx]
            )
            winch_force[winch.idx] ~ F
        ]
        defaults = [
            defaults
            tether_length[winch.idx] => winch.tether_length
            tether_vel[winch.idx] => 0
        ]
    end
    return eqs, defaults, guesses, set_values, tether_kite_force, tether_kite_torque_b
end

function create_sys!(s::KPSQ, system::PointMassSystem, wing::KiteWing; I_p, R_b_p, Q_p_b, init_Q_p_w, init_kite_pos, init=false)
    eqs = []
    defaults = Pair{Num, Real}[]
    guesses = Pair{Num, Real}[]
    tether_kite_force = zeros(Num, 3)
    tether_kite_torque_b = zeros(Num, 3)

    @parameters begin
        measured_wind_dir_gnd = s.measure.wind_dir_gnd
        # measured_sphere_pos[1:2, 1:2] = s.measure.sphere_pos
        # measured_sphere_vel[1:2, 1:2] = s.measure.sphere_vel
        # measured_sphere_acc[1:2, 1:2] = s.measure.sphere_acc
        # measured_tether_length[1:3] = s.measure.tether_length
        # measured_tether_vel[1:3]    = s.measure.tether_vel
        # measured_tether_acc[1:3]    = s.measure.tether_acc
    end
    @variables begin
        # potential differential variables
        kite_pos(t)[1:3] # xyz pos of kite in world frame
        kite_vel(t)[1:3]
        kite_acc(t)[1:3]
        distance(t)
        distance_vel(t)
        distance_acc(t)
        ω_p(t)[1:3] # turn rate in principal frame
        ω_b(t)[1:3] # turn rate in body frame
        α_p(t)[1:3] # angular acceleration in principal frame
        α_b(t)[1:3] # angular acceleration in body frame

        # rotations and frames
        R_b_w(t)[1:3, 1:3] # rotation of the kite body frame relative to the world frame
        R_p_w(t)[1:3, 1:3] # rotation of the kite principal frame relative to the world frame
        Q_p_w(t)[1:4] # quaternion orientation of the kite principal frame relative to the world frame
        Q_b_w(t)[1:4] # quaternion orientation of the kite body frame relative to the world frame
        Q_vel(t)[1:4] # quaternion rate of change
        e_x(t)[1:3]
        e_y(t)[1:3]
        e_z(t)[1:3]

        # rest: forces, torques, vectors and scalar values
        torque_p(t)[1:3] # torque in principal frame
        torque_b(t)[1:3] # torque in body frame
        total_kite_force(t)[1:3]
        aero_kite_force(t)[1:3]
        rho_kite(t)
        wind_vec_gnd(t)[1:3]
        wind_vel_kite(t)[1:3]
        va_kite(t)[1:3]
        va_kite_b(t)[1:3]
    end

    q_inf = (0.5 * rho_kite * norm(va_kite_b)^2)

    function diff_eqs!()
        @parameters begin
            force_coefficients[1:3] = s.vsm_solver.sol.force_coefficients
            torque_coefficients[1:3] = s.vsm_solver.sol.moment_coefficients
        end
        Q_b_p = quaternion_conjugate(Q_p_b)
        Ω = [0       -ω_p[1]  -ω_p[2]  -ω_p[3];
            ω_p[1]    0        ω_p[3]  -ω_p[2];
            ω_p[2]   -ω_p[3]   0        ω_p[1];
            ω_p[3]    ω_p[2]  -ω_p[1]   0]
        eqs = [
            eqs
            [D(Q_p_w[i]) ~ Q_vel[i] for i in 1:4]
            [Q_vel[i] ~ 0.5 * sum(Ω[i, j] * Q_p_w[j] for j in 1:4) for i in 1:4]
            [R_b_w[:, i] ~ quaternion_to_rotation_matrix(quaternion_multiply(Q_p_w, Q_b_p))[:, i] for i in 1:3] # https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation#Performance_comparisons
            [R_p_w[:, i] ~ quaternion_to_rotation_matrix(Q_p_w)[:, i] for i in 1:3]
            D(ω_p) ~ α_p
            ω_b ~ R_b_p' * ω_p
            α_p[1] ~ (torque_p[1] + (I_p[2] - I_p[3]) * ω_p[2] * ω_p[3]) / I_p[1]
            α_p[2] ~ (torque_p[2] + (I_p[3] - I_p[1]) * ω_p[3] * ω_p[1]) / I_p[2]
            α_p[3] ~ (torque_p[3] + (I_p[1] - I_p[2]) * ω_p[1] * ω_p[2]) / I_p[3]
            α_b ~ R_b_p' * α_p
            torque_b ~ torque_coefficients * q_inf * s.aero.projected_area + tether_kite_torque_b
            torque_p ~ R_b_p * torque_b
            
            [D(kite_pos[i]) ~ kite_vel[i] for i in 1:3]
            [D(kite_vel[i]) ~ kite_acc[i] for i in 1:3]
            aero_kite_force ~ (R_b_w * force_coefficients) * q_inf * s.aero.projected_area
            kite_acc        ~ (tether_kite_force + aero_kite_force) / s.set.mass

            distance            ~ norm(kite_pos)
            distance_vel        ~ kite_vel ⋅ normalize(kite_pos)
            distance_acc        ~ kite_acc ⋅ normalize(kite_pos)    
        ]
        defaults = [
            defaults
            [Q_p_w[i] => init_Q_p_w[i] for i in 1:4]
            [ω_p[i] => 0 for i in 1:3]
            [kite_pos[i] => init_kite_pos[i] for i in 1:3]
            [kite_vel[i] => 0 for i in 1:3]
        ]
        return nothing
    end
    
    function scalar_eqs!()
        @parameters wind_scale_gnd = s.set.v_wind
        eqs = [
            eqs
            e_x     ~ R_b_w * [1, 0, 0]
            e_y     ~ R_b_w * [0, 1, 0]
            e_z     ~ R_b_w * [0, 0, 1]
            rho_kite        ~ calc_rho(s.am, kite_pos[3])
            wind_vec_gnd ~ wind_scale_gnd * rotate_around_z([1, 0, 0], measured_wind_dir_gnd)
            wind_vel_kite  ~ AtmosphericModels.calc_wind_factor(s.am, kite_pos[3], s.set.profile_law) * wind_vec_gnd
            va_kite ~ wind_vel_kite - kite_vel
            va_kite_b ~ R_b_w' * va_kite
        ]
        @variables begin
            heading_y(t)
            azimuth(t)
            azimuth_vel(t)
            azimuth_acc(t)
            elevation(t)
            elevation_vel(t)
            elevation_acc(t)
            x_acc(t)
            y_acc(t)
            sphere_pos(t)[1:2, 1:2] # TODO: add equations for these, and think of a good measurement system, maybe rolling window?
            sphere_vel(t)[1:2, 1:2]
            sphere_acc(t)[1:2, 1:2]
        end

        x, y, z = kite_pos
        x´, y´, z´ = kite_vel
        x´´, y´´, z´´ = kite_acc

        eqs = [
            eqs
            heading_y       ~ calc_heading_y(-e_x)

            elevation           ~ atan(z / x)
            # elevation_vel = d/dt(atan(z/x)) = (x*ż' - z*ẋ')/(x^2 + z^2) according to wolframalpha
            elevation_vel       ~ (x*z´ - z*x´) / 
                                    (x^2 + z^2)
            elevation_acc       ~ ((x^2 + z^2)*(x*z´´ - z*x´´) + 2(z*x´ - x*z´)*(x*x´ + z*z´))/(x^2 + z^2)^2
            azimuth             ~ atan(y / x)
            # azimuth_vel = d/dt(atan(y/x)) = (-y*x´ + x*y´)/(x^2 + y^2) # TODO: check if correct
            azimuth_vel         ~ (-y*x´ + x*y´) / 
                                    (x^2 + y^2)
            azimuth_acc         ~ ((x^2 + y^2)*(-y*x´´ + x*y´´) + 2(y*x´ - x*y´)*(x*x´ + y*y´))/(x^2 + y^2)^2
            x_acc               ~ kite_acc ⋅ e_x
            y_acc               ~ kite_acc ⋅ e_y
        ]
        return nothing
    end

    eqs, defaults, guesses, set_values, tether_kite_force, tether_kite_torque_b = 
        force_eqs!(s, system, eqs, defaults, guesses; tether_kite_force, tether_kite_torque_b, R_b_w, kite_pos, q_inf, wind_vec_gnd)
    diff_eqs!()
    scalar_eqs!()
    
    # te_I = (1/3 * (s.set.mass/8) * te_length^2)
    # # -damping / I * ω = α_damping
    # # solve for c: (c * (k*m/s^2) / (k*m^2)) * (m/s)=m/s^2 in wolframalpha
    # # damping should be in N*m*s
    # rot_damping = 0.1s.damping * te_length

    # eqs = [
    #     eqs
    #     trailing_edge_α[1] ~ (force[:, s.i_A]) ⋅ e_te_A * te_length / te_I - (rot_damping[1] / te_I) * trailing_edge_ω[1]
    #     trailing_edge_α[2] ~ (force[:, s.i_B]) ⋅ e_te_B * te_length / te_I - (rot_damping[2] / te_I) * trailing_edge_ω[2]
    # ]
    
    eqs = Symbolics.scalarize.(reduce(vcat, Symbolics.scalarize.(eqs)))

    # if !init
    #     discrete_events = [true => [Q_p_w[i] ~ normalize(Q_p_w)[i] for i in 1:4]]
    # else
    #     discrete_events = [true => [Q_b_w[i] ~ normalize(Q_b_w)[i] for i in 1:4]]
    # end
    
    # @named sys = ODESystem(eqs, t; discrete_events)
    @info "Creating ODESystem"
    @time @named sys = ODESystem(eqs, t)
    return sys, collect(set_values), defaults, guesses
end

function model!(s::KPSQ; init=false)
    I_p, I_b, R_b_p, Q_p_b = calc_inertia(s.wing)
    init_Q_p_w, R_b_w = measure_to_q(s.measure, R_b_p)
    @show R_b_w

    s.point_system = PointMassSystem(s, s.wing)
    init_kite_pos = init!(s.point_system, s, R_b_w)
    @show init_kite_pos

    VortexStepMethod.set_va!(s.aero, [s.set.v_wind, 0., 0.])
    VortexStepMethod.solve!(s.vsm_solver, s.aero)
    
    sys, inputs, defaults, guesses = create_sys!(s, s.point_system, s.wing; I_p, R_b_p, Q_p_b, init_Q_p_w, init_kite_pos, init)
    @info "Simplifying the system"
    @time sys = structural_simplify(sys; simplify=false)
    return sys, defaults, guesses
end
