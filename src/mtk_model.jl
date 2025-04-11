# ==================== mtk model functions ================================================
# Implementation of the ram air kite model using ModelingToolkit.jl

function calc_speed_acc(winch::AsyncMachine, tether_vel, norm_, set_speed)
    calc_acceleration(winch, tether_vel, norm_; set_speed, set_torque=nothing, use_brake=false) # TODO: add brake setting
end
function calc_moment_acc(winch::TorqueControlledMachine, tether_vel, norm_, set_torque)
    calc_acceleration(winch, tether_vel, norm_; set_speed=nothing, set_torque, use_brake=false)
end
@register_symbolic calc_speed_acc(winch::AsyncMachine, tether_vel, norm_, set_speed)
@register_symbolic calc_moment_acc(winch::TorqueControlledMachine, tether_vel, norm_, set_torque)

function sym_interp(interp::Function, aoa, trailing_edge_angle)
    return interp(rad2deg(aoa), rad2deg(trailing_edge_angle-aoa)) # TODO: register callable struct https://docs.sciml.ai/Symbolics/dev/manual/functions/#Symbolics.@register_array_symbolic
end
@register_symbolic sym_interp(interp::Function, aoa, trailing_edge_angle)

function sym_normalize(vec)
    return vec / norm(vec)
end
@register_symbolic sym_normalize(vec)

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
        R_b_w, kite_pos, kite_vel, wind_vec_gnd, group_aero_moment, twist_angle, steady)

    @parameters acc_multiplier = 1

    points, groups, segments, pulleys, tethers, winches = 
        system.points, system.groups, system.segments, system.pulleys, system.tethers, system.winches
    
    # ==================== POINTS ==================== #
    @variables begin
        pos(t)[1:3, eachindex(points)]
        vel(t)[1:3, eachindex(points)]
        acc(t)[1:3, eachindex(points)]
        point_force(t)[1:3, eachindex(points)]
        tether_kite_force(t)[1:3, eachindex(points)]
        tether_kite_moment(t)[1:3, eachindex(points)]
        tether_r(t)[1:3, eachindex(points)]
        point_mass(t)[eachindex(points)]

        spring_force_vec(t)[1:3, eachindex(segments)]
        drag_force(t)[1:3, eachindex(segments)]
        l0(t)[eachindex(segments)]
    end
    for point in points
        F::Vector{Num} = zeros(Num, 3)
        mass = 0.0
        in_bridle = false
        for segment in segments
            if point.idx in segment.points
                mass_per_meter = s.set.rho_tether * π * (segment.diameter/2)^2
                inverted = segment.points[2] == point.idx
                if inverted
                    F .-= spring_force_vec[:, segment.idx]
                else
                    F .+= spring_force_vec[:, segment.idx]
                end
                mass += mass_per_meter * l0[segment.idx] / 2
                F .+= 0.5drag_force[:, segment.idx]

                if segment.type == BRIDLE
                    in_bridle = true
                end
            end
        end
        eqs = [
            eqs
            point_mass[point.idx] ~ mass
            point_force[:, point.idx]  ~ F
        ]

        winch_point = false
        for tether in tethers
            if point.idx == tether.winch_point
                winch_point = true
                break
            end
        end

        if point.type !== KITE
            eqs = [
                eqs
                tether_kite_force[:, point.idx] ~ zeros(3)
                tether_kite_moment[:, point.idx] ~ zeros(3)
            ]
        end

        if point.type == KITE
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

            group = groups[group_idx]
            if point.idx == group.points[group.fixed_index]
                pos_b = point.pos_b
            else
                fixed_pos = points[group.points[group.fixed_index]].pos_b
                chord_b = point.pos_b - fixed_pos
                normal = chord_b × group.y_airf
                pos_b = fixed_pos + cos(twist_angle[group_idx]) * chord_b - sin(twist_angle[group_idx]) * normal
            end
            eqs = [
                eqs
                tether_kite_force[:, point.idx] ~ point_force[:, point.idx]
                tether_r[:, point.idx] ~ pos[:, point.idx] - kite_pos
                tether_kite_moment[:, point.idx] ~ tether_r[:, point.idx] × point_force[:, point.idx]
                pos[:, point.idx]    ~ kite_pos + R_b_w * pos_b
                vel[:, point.idx]    ~ zeros(3)
                acc[:, point.idx]    ~ zeros(3)
            ]
        elseif point.type == WINCH
            eqs = [
                eqs
                pos[:, point.idx]    ~ point.pos_w
                vel[:, point.idx]    ~ zeros(3)
                acc[:, point.idx]    ~ zeros(3)
            ]
        elseif point.type == DYNAMIC
            p = pos[:, point.idx]
            n = sym_normalize(kite_pos)
            n = n * (p ⋅ n)
            r = (p - n) # https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line#Vector_formulation
            @parameters bridle_damp = 100
            @parameters measured_ω_z = 0.6
            eqs = [
                eqs
                D(pos[:, point.idx]) ~ vel[:, point.idx]
                D(vel[:, point.idx]) ~ acc_multiplier * acc[:, point.idx] - in_bridle * bridle_damp * (vel[:, point.idx] - kite_vel)
                acc[:, point.idx]    ~ point_force[:, point.idx] / mass + 
                                        [0, 0, -G_EARTH] + 
                                        ifelse.(steady==true, r * norm(measured_ω_z)^2, zeros(3)) # TODO: add other steady accelerations
            ]
            defaults = [
                defaults
                [pos[j, point.idx] => point.pos_w[j] for j in 1:3]
                [vel[j, point.idx] => 0 for j in 1:3]
            ]
        elseif point.type == STATIC
            eqs = [
                eqs
                vel[:, point.idx]    ~ zeros(3)
                acc[:, point.idx]    ~ zeros(3)
                acc[:, point.idx]    ~ point_force[:, point.idx] / mass + [0, 0, -G_EARTH]
            ]
            guesses = [
                guesses
                [acc[j, point.idx] => 0 for j in 1:3]
                [pos[j, point.idx] => point.pos_w[j] for j in 1:3]
                [point_force[j, point.idx] => 0 for j in 1:3]
            ]
        else
            throw(ArgumentError("Unknown point type: $(typeof(point))"))
        end
    end

    # ==================== GROUPS ==================== #
    @variables begin
        trailing_edge_angle(t)[eachindex(groups)] # angle left / right
        trailing_edge_ω(t)[eachindex(groups)] # angular rate
        trailing_edge_α(t)[eachindex(groups)] # angular acc
        free_twist_angle(t)[eachindex(groups)]
        twist_ω(t)[eachindex(groups)] # angular rate
        twist_α(t)[eachindex(groups)] # angular acc
        group_tether_moment(t)[eachindex(groups)]
        tether_moment(t)[eachindex(groups), 1:length(groups[1].points)-1]
    end
    
    for group in groups

        x_airf = normalize(group.chord)
        init_z_airf = x_airf × group.y_airf
        z_airf = x_airf * sin(twist_angle[group.idx]) + init_z_airf * cos(twist_angle[group.idx])
        moving_points = filter(p -> p != group.points[group.fixed_index], group.points)
        for (i, point_idx) in enumerate(moving_points)
            r = (points[point_idx].pos_b - points[group.points[group.fixed_index]].pos_b) ⋅ normalize(group.chord)
            eqs = [
                eqs
                tether_moment[group.idx, i] ~ r * (point_force[:, point_idx] ⋅ (R_b_w * -z_airf))
            ]
        end
        
        inertia = 1/3 * (s.set.mass/length(groups)) * (norm(group.chord))^2 # plate inertia around leading edge
        @assert !(inertia ≈ 0.0)
        @parameters twist_damp = s.set.quasi_static ? 200 : 100
        eqs = [
            eqs
            group_tether_moment[group.idx] ~ sum(tether_moment[group.idx, :])
            twist_α[group.idx] ~ (group_aero_moment[group.idx] + group_tether_moment[group.idx]) / inertia
            twist_angle[group.idx] ~ clamp(free_twist_angle[group.idx], -π/2, π/2)
        ]
        if group.type == DYNAMIC
            eqs = [
                eqs
                D(free_twist_angle[group.idx]) ~ twist_ω[group.idx]
                D(twist_ω[group.idx]) ~ twist_α[group.idx] - twist_damp * twist_ω[group.idx]
            ]
            defaults = [
                defaults
                free_twist_angle[group.idx] => 0
                twist_ω[group.idx] => 0
            ]
        elseif group.type == STATIC
            eqs = [
                eqs
                twist_ω[group.idx] ~ 0
                twist_α[group.idx] ~ 0
            ]
            guesses = [
                guesses
                free_twist_angle[group.idx] => 0
                twist_angle[group.idx] => 0
            ]
        else
            throw(ArgumentError("Wrong group type."))
        end
    end

    # ==================== SEGMENTS ==================== #
    @variables begin
        segment_vec(t)[1:3, eachindex(segments)]
        unit_vector(t)[1:3, eachindex(segments)]
        len(t)[eachindex(segments)]
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

        tether_length(t)[eachindex(winches)]
    end
    for segment in segments
        p1, p2 = segment.points[1], segment.points[2]
        guesses = [
            guesses
            [segment_vec[i, segment.idx] => points[p2].pos_w[i] - points[p1].pos_w[i] for i in 1:3]
        ]

        if segment.type == BRIDLE
            in_pulley = 0
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
        elseif segment.type == POWER || segment.type == STEERING
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

        stiffness_m = s.set.e_tether * (segment.diameter/2)^2 * pi
        (segment.type == BRIDLE) && (compression_frac = 0.1)
        (segment.type == POWER) && (compression_frac = 0.1)
        (segment.type == STEERING) && (compression_frac = 0.1)
        
        @parameters stiffness_frac = 0.1
        (segment.type == BRIDLE) && (stiffness_m = stiffness_frac * stiffness_m)

        damping_m = (s.set.damping / s.set.c_spring) * stiffness_m
        
        eqs = [
            eqs
            # spring force equations
            segment_vec[:, segment.idx]  ~ pos[:, p2] - pos[:, p1]
            len[segment.idx]             ~ norm(segment_vec[:, segment.idx])
            unit_vector[:, segment.idx]  ~ segment_vec[:, segment.idx]/len[segment.idx]
            rel_vel[:, segment.idx]      ~ vel[:, p1] - vel[:, p2]
            spring_vel[segment.idx]      ~ rel_vel[:, segment.idx] ⋅ unit_vector[:, segment.idx]
            stiffness[segment.idx]       ~ ifelse(len[segment.idx] > l0[segment.idx],
                                        stiffness_m / len[segment.idx],
                                        compression_frac * stiffness_m / len[segment.idx]
            )
            damping[segment.idx]         ~ damping_m / len[segment.idx]
            spring_force[segment.idx] ~  (stiffness[segment.idx] * (len[segment.idx] - l0[segment.idx]) - 
                            damping[segment.idx] * spring_vel[segment.idx])
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
    @parameters pulley_damp = 20
    for pulley in pulleys
        segment = segments[pulley.segments[1]]
        mass_per_meter = s.set.rho_tether * π * (segment.diameter/2)^2
        mass = pulley.sum_length * mass_per_meter
        @assert !(mass ≈ 0.0)
        eqs = [
            eqs
            pulley_force[pulley.idx]    ~ spring_force[pulley.segments[1]] - spring_force[pulley.segments[2]]
            pulley_acc[pulley.idx]      ~ pulley_force[pulley.idx] / mass
        ]
        if pulley.type == DYNAMIC
            eqs = [
                eqs
                D(pulley_l0[pulley.idx])  ~ pulley_vel[pulley.idx]
                D(pulley_vel[pulley.idx]) ~ acc_multiplier * pulley_acc[pulley.idx] - pulley_damp * pulley_vel[pulley.idx]
            ]
            defaults = [
                defaults
                pulley_l0[pulley.idx] => segments[pulley.segments[1]].l0
                pulley_vel[pulley.idx] => 0
            ]
        elseif pulley.type == STATIC
            eqs = [
                eqs 
                pulley_vel[pulley.idx] ~ 0
                pulley_acc[pulley.idx] ~ 0
            ]
            guesses = [
                guesses
                pulley_l0[pulley.idx] => segments[pulley.segments[1]].l0
            ]
        else
            throw(ArgumentError("Wrong pulley type"))
        end
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
            D(tether_vel[winch.idx]) ~ ifelse(steady==true, 0, tether_acc[winch.idx])

            tether_acc[winch.idx] ~ calc_moment_acc( # TODO: moment and speed control
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
    return eqs, defaults, guesses, tether_kite_force, tether_kite_moment
end

function diff_eqs!(s, eqs, defaults; tether_kite_force, tether_kite_moment, aero_force_b, 
    aero_moment_b, ω_b, R_b_w, kite_pos, kite_vel, kite_acc, init_Q_b_w, init_kite_pos, steady
)
    @variables begin
        # potential differential variables
        kite_acc_b(t)[1:3]
        α_b(t)[1:3] # angular acceleration in principal frame

        # rotations and frames
        Q_b_w(t)[1:4] # quaternion orientation of the kite body frame relative to the world frame
        Q_vel(t)[1:4] # quaternion rate of change

        # rest: forces, moments, vectors and scalar values
        moment_b(t)[1:3] # moment in principal frame
        total_tether_kite_force(t)[1:3]
        total_tether_kite_moment(t)[1:3]
    end
    @parameters ω_damp = 150

    Ω = [0       -ω_b[1]  -ω_b[2]  -ω_b[3];
        ω_b[1]    0        ω_b[3]  -ω_b[2];
        ω_b[2]   -ω_b[3]   0        ω_b[1];
        ω_b[3]    ω_b[2]  -ω_b[1]   0]

    I_b = [s.wing.inertia_tensor[1,1], s.wing.inertia_tensor[2,2], s.wing.inertia_tensor[3,3]]
    @assert !(s.set.mass ≈ 0)
    eqs = [
        eqs
        [D(Q_b_w[i]) ~ Q_vel[i] for i in 1:4]
        [Q_vel[i] ~ 0.5 * sum(Ω[i, j] * Q_b_w[j] for j in 1:4) for i in 1:4]
        [R_b_w[:, i] ~ quaternion_to_rotation_matrix(Q_b_w)[:, i] for i in 1:3]

        D(ω_b[1]) ~ α_b[1]
        D(ω_b[2]) ~ α_b[2] - ω_damp * ω_b[2]
        D(ω_b[3]) ~ ifelse(steady==true, 0, α_b[3])

        α_b[1] ~ (moment_b[1] + (I_b[2] - I_b[3]) * ω_b[2] * ω_b[3]) / I_b[1]
        α_b[2] ~ (moment_b[2] + (I_b[3] - I_b[1]) * ω_b[3] * ω_b[1]) / I_b[2]
        α_b[3] ~ (moment_b[3] + (I_b[1] - I_b[2]) * ω_b[1] * ω_b[2]) / I_b[3]
        total_tether_kite_moment ~ [sum(tether_kite_moment[i, :]) for i in 1:3]
        moment_b ~ aero_moment_b + R_b_w' * total_tether_kite_moment
        
        D(kite_pos) ~ kite_vel
        D(kite_vel) ~ kite_acc
        kite_acc ~ R_b_w * [ifelse(steady==true, 0, kite_acc_b[1]),
                            ifelse(steady==true, 0, kite_acc_b[2]),
                            kite_acc_b[3]]
        kite_acc_b        ~ (R_b_w' * total_tether_kite_force + aero_force_b) / s.set.mass
        total_tether_kite_force ~ [sum(tether_kite_force[i, :]) for i in 1:3]
    ]
    defaults = [
        defaults
        [Q_b_w[i] => init_Q_b_w[i] for i in 1:4]
        [ω_b[i] => 0 for i in 1:3]
        [kite_pos[i] => init_kite_pos[i] for i in 1:3]
        [kite_vel[i] => 0 for i in 1:3]
    ]
    return eqs, defaults
end


function scalar_eqs!(s, eqs, measure; R_b_w, wind_vec_gnd, va_kite_b, kite_pos, kite_vel, kite_acc)
    @parameters wind_scale_gnd = s.set.v_wind
    @parameters measured_wind_dir_gnd = measure.wind_dir_gnd
    @variables begin
        e_x(t)[1:3]
        e_y(t)[1:3]
        e_z(t)[1:3]
        wind_vel_kite(t)[1:3]
        va_kite(t)[1:3]
    end
    eqs = [
        eqs
        e_x     ~ R_b_w * [1, 0, 0]
        e_y     ~ R_b_w * [0, 1, 0]
        e_z     ~ R_b_w * [0, 0, 1]
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
        course(t)
        x_acc(t)
        y_acc(t)
        sphere_pos(t)[1:2, 1:2]
        sphere_vel(t)[1:2, 1:2]
        sphere_acc(t)[1:2, 1:2]
    end

    x, y, z = kite_pos
    x´, y´, z´ = kite_vel
    x´´, y´´, z´´ = kite_acc

    eqs = [
        eqs
        heading_y       ~ atan(e_x[2]/e_x[1])

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
        course              ~ atan(-azimuth_vel, elevation_vel)
    ]
    return eqs
end

function linear_vsm_eqs!(s, eqs; aero_force_b, aero_moment_b, group_aero_moment, init_va, twist_angle, va_kite_b, ω_b)
    sol = s.vsm_solver.sol
    @assert length(s.point_system.groups) == length(sol.group_moment_dist)

    y_ = [zeros(4); init_va; zeros(3)]
    jac_, x_ = VortexStepMethod.linearize(
        s.vsm_solver, 
        s.aero, 
        y_;
        theta_idxs=1:4, 
        va_idxs=5:7, 
        omega_idxs=8:10,
        moment_frac=s.bridle_fracs[s.point_system.groups[1].fixed_index])

    @parameters begin
        last_y[eachindex(y_)] = y_
        last_x[eachindex(x_)] = x_
        vsm_jac[eachindex(x_), eachindex(y_)] = jac_
    end

    @variables begin
        y(t)[eachindex(y_)]
        dy(t)[eachindex(y_)]
    end

    eqs = [
        eqs
        y ~ [twist_angle; va_kite_b; ω_b]
        dy ~ y - last_y
        [aero_force_b; aero_moment_b; group_aero_moment] ~ last_x + vsm_jac * dy
    ]

    return eqs
end

function init_unknowns_vec!(
    s::RamAirKite, 
    system::PointMassSystem, 
    vec::Vector{SimFloat},
    init_Q_b_w,
    init_kite_pos;
    non_observed=true
)
    !s.set.quasi_static && non_observed && (length(vec) != length(s.integrator.u)) && 
        throw(ArgumentError("Unknowns of length $(length(s.integrator.u)) but vector provided of length $(length(vec))"))
        
    @unpack points, groups, segments, pulleys, winches = system
    vec_idx = 1
    
    if non_observed
        for point in points
            if point.type == DYNAMIC
                for i in 1:3
                    vec[vec_idx] = point.pos_w[i]
                    vec_idx += 1
                end
                for i in 1:3 # TODO: add speed to init
                    vec[vec_idx] = 0.0
                    vec_idx += 1
                end
            end
        end

        for group in groups
            if group.type == DYNAMIC
                vec[vec_idx] = 0
                vec_idx += 1
                vec[vec_idx] = 0
                vec_idx += 1
            end
        end

        for pulley in pulleys
            if pulley.type == DYNAMIC
                vec[vec_idx] = segments[pulley.segments[1]].l0
                vec_idx += 1
                vec[vec_idx] = 0
                vec_idx += 1
            end
        end
    end

    for winch in winches
        vec[vec_idx] = winch.tether_length
        vec_idx += 1
        vec[vec_idx] = 0
        vec_idx += 1
    end

    for i in 1:4
        vec[vec_idx] = init_Q_b_w[i]
        vec_idx += 1
    end
    for i in 1:3
        vec[vec_idx] = 0
        vec_idx += 1
    end
    for i in 1:3
        vec[vec_idx] = init_kite_pos[i]
        vec_idx += 1
    end
    for i in 1:3
        vec[vec_idx] = 0
        vec_idx += 1
    end
    non_observed && (vec_idx-1 != length(vec)) && 
        throw(ArgumentError("Unknowns vec is of length $(length(vec)) but the last index is $(vec_idx-1)"))
    nothing
end

function get_unknowns(s::RamAirKite)
    vec = Num[]
    vec = get_stiff_unknowns(s, vec)
    vec = get_nonstiff_unknowns(s, vec)
    !s.set.quasi_static && (length(vec) != length(s.integrator.u)) &&
        throw(ArgumentError("Integrator unknowns of length $(length(s.integrator.u)) should equal vec of length $(length(vec))"))
    return vec
end

function get_stiff_unknowns(s, vec=Num[])
    @unpack points, groups, segments, pulleys, winches = s.point_system
    for point in points
        for i in 1:3
            point.type == DYNAMIC && push!(vec, s.sys.pos[i, point.idx])
        end
        for i in 1:3 # TODO: add speed to init
            point.type == DYNAMIC && push!(vec, s.sys.vel[i, point.idx])
        end
    end
    for group in groups
        group.type == DYNAMIC && push!(vec, s.sys.free_twist_angle[group.idx])
        group.type == DYNAMIC && push!(vec, s.sys.twist_ω[group.idx])
    end
    for pulley in pulleys
        pulley.type == DYNAMIC && push!(vec, s.sys.pulley_l0[pulley.idx])
        pulley.type == DYNAMIC && push!(vec, s.sys.pulley_vel[pulley.idx])
    end
    return vec
end

function get_nonstiff_unknowns(s, vec=Num[])
    @unpack points, groups, segments, pulleys, winches = s.point_system
    for winch in winches
        push!(vec, s.sys.tether_length[winch.idx])
        push!(vec, s.sys.tether_vel[winch.idx])
    end
    [push!(vec, s.sys.Q_b_w[i]) for i in 1:4]
    [push!(vec, s.sys.ω_b[i]) for i in 1:3]
    [push!(vec, s.sys.kite_pos[i]) for i in 1:3]
    [push!(vec, s.sys.kite_vel[i]) for i in 1:3]
    return vec
end

function create_sys!(s::RamAirKite, system::PointMassSystem, measure::Measurement; init_Q_b_w, init_kite_pos, init_va)
    eqs = []
    defaults = Pair{Num, Real}[]
    guesses = Pair{Num, Real}[]

    @parameters begin
        steady = false
        # measured_sphere_pos[1:2, 1:2] = measure.sphere_pos
        # measured_sphere_vel[1:2, 1:2] = measure.sphere_vel
        # measured_sphere_acc[1:2, 1:2] = measure.sphere_acc
        # measured_tether_length[1:3] = measure.tether_length
        # measured_tether_vel[1:3]    = measure.tether_vel
        # measured_tether_acc[1:3]    = measure.tether_acc
    end
    @variables begin
        # potential differential variables
        kite_pos(t)[1:3] # xyz pos of kite in world frame
        kite_vel(t)[1:3]
        kite_acc(t)[1:3]
        ω_b(t)[1:3] # turn rate in principal frame

        # rotations and frames
        R_b_w(t)[1:3, 1:3] # rotation of the kite body frame relative to the world frame

        # rest: forces, moments, vectors and scalar values
        aero_force_b(t)[1:3]
        aero_moment_b(t)[1:3]
        twist_angle(t)[eachindex(system.groups)]
        group_aero_moment(t)[eachindex(system.groups)]
        wind_vec_gnd(t)[1:3]
        va_kite_b(t)[1:3]
    end

    eqs, defaults, guesses, tether_kite_force, tether_kite_moment = 
        force_eqs!(s, system, eqs, defaults, guesses; 
            R_b_w, kite_pos, kite_vel, wind_vec_gnd, group_aero_moment, twist_angle, steady)
    eqs = linear_vsm_eqs!(s, eqs; aero_force_b, aero_moment_b, group_aero_moment, init_va, twist_angle, va_kite_b, ω_b)
    eqs, defaults = diff_eqs!(s, eqs, defaults; tether_kite_force, tether_kite_moment, aero_force_b, aero_moment_b, 
        ω_b, R_b_w, kite_pos, kite_vel, kite_acc, init_Q_b_w, init_kite_pos, steady)
    eqs = scalar_eqs!(s, eqs, measure; R_b_w, wind_vec_gnd, va_kite_b, kite_pos, kite_vel, kite_acc)
    
    # te_I = (1/3 * (s.set.mass/8) * te_length^2)
    # # -damping / I * ω = α_damping
    # # solve for c: (c * (k*m/s^2) / (k*m^2)) * (m/s)=m/s^2 in wolframalpha
    # # damping should be in N*m*s
    # rot_damping = 0.1s.damping * te_length

    # eqs = [
    #     eqs
    #     trailing_edge_α[1] ~ (force[:, s.i_A]) ⋅ e_te_A * te_length / te_I - (rot_damping[1] / te_I) * trailing_edge_ω[1] # TODO: add trailing edge
    #     trailing_edge_α[2] ~ (force[:, s.i_B]) ⋅ e_te_B * te_length / te_I - (rot_damping[2] / te_I) * trailing_edge_ω[2]
    # ]
    
    eqs = Symbolics.scalarize.(reduce(vcat, Symbolics.scalarize.(eqs)))

    # discrete_events = [
    #     true => [
    #         [Q_b_w[i] ~ normalize(Q_b_w)[i] for i in 1:4]
    #         [twist_angle[i] ~ clamp(twist_angle[i], -π/2, π/2) for i in eachindex(s.point_system.groups)]
    #         ]
    #     ]

    @info "Creating ODESystem"
    # @named sys = ODESystem(eqs, t; discrete_events)
    @time @named sys = ODESystem(eqs, t)
    return sys, defaults, guesses
end

