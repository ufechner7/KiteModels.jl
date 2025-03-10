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

function mean(vec)
    return sum(vec) / length(vec)
end

function convert_pos_vel(s::KPSQ, pos_, vel_)
    pos = Array{Union{Nothing, SimFloat}}(nothing, 3, s.i_A)
    vel = Array{Union{Nothing, SimFloat}}(nothing, 3, s.i_A)
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


function create_point_mass_system!(s::KPSQ, wing::KiteWing)
    # TODO: move as much of the code as possible from create_point_mass_system to other places, to make model creation easier.
    # 1. move bridle gamma calculation
    # 2. ...

    points = AbstractPoint[]
    groups = KitePointGroup[]
    segments = Segment[]
    pulleys = Pulley[]
    tethers = Tether[]

    attach_points = AbstractPoint[]
    
    bridle_gammas, bridle_limits = find_bridle_gammas!(s, wing)

    function create_bridle(gammas, limits)
        i_pnt = length(points) # last point idx
        i_seg = length(segments) # last segment idx
        i_pul = length(pulleys) # last pulley idx

        i = 1
        for (gamma, limit) in zip(gammas, limits) # 2 gammas with 2 pairs of limits
            le_pos = [wing.le_interp[i](gamma) for i in 1:3]
            chord = [wing.te_interp[i](gamma) for i in 1:3] .- le_pos
            y_panel = normalize(le_pos .- [wing.le_interp[i](gamma+0.01) for i in 1:3])
            fixed_pos = le_pos .+ chord .* s.bridle_fracs[2]
            point_idxs = Int16[]
            for frac in s.bridle_fracs # 4 fracs
                pos = le_pos .+ chord .* frac
                points = [points; KitePoint(i+i_pnt, pos, fixed_pos, y_panel)]
                push!(point_idxs, points[end].idx)
                i += 1
            end
            
            i_grp = 1 + length(groups)
            y_lim = (wing.le_interp[2](limit[1]), wing.le_interp[2](limit[2]))
            @show y_lim
            groups = [groups; KitePointGroup(i_grp, point_idxs, y_lim)]
        end

        mean_le = [wing.le_interp[i](mean(gammas)) for i in 1:3]
        chord_length = norm([wing.te_interp[i](mean(gammas)) for i in 1:3] .- mean_le)
        xs = s.bridle_fracs .* chord_length
        bridle_top = mean_le .+ [0, 0, -3]

        points = [
            points
            Point(9+i_pnt, bridle_top .+ [xs[1], 0, 0])
            Point(10+i_pnt, bridle_top .+ [xs[2], 0, 0])
            Point(11+i_pnt, bridle_top .+ [xs[3], 0, 0])
            Point(12+i_pnt, bridle_top .+ [xs[4], 0, 0])

            Point(13+i_pnt, bridle_top .+ [xs[2], 0, -1])

            Point(14+i_pnt, bridle_top .+ [xs[1], 0, -2])
            Point(15+i_pnt, bridle_top .+ [xs[3], 0, -2])

            Point(16+i_pnt, bridle_top .+ [xs[1], 0, -5])
            Point(17+i_pnt, bridle_top .+ [xs[3], 0, -5])
        ]
        segments = [
            segments
            Segment(1+i_seg, (1+i_pnt, 9+i_pnt), 2, BRIDLE)
            Segment(2+i_seg, (2+i_pnt, 10+i_pnt), 2, BRIDLE)
            Segment(3+i_seg, (3+i_pnt, 11+i_pnt), 2, BRIDLE)
            Segment(4+i_seg, (4+i_pnt, 12+i_pnt), 2, BRIDLE)

            Segment(5+i_seg, (5+i_pnt, 9+i_pnt), 2, BRIDLE)
            Segment(6+i_seg, (6+i_pnt, 10+i_pnt), 2, BRIDLE)
            Segment(7+i_seg, (7+i_pnt, 11+i_pnt), 2, BRIDLE)
            Segment(8+i_seg, (8+i_pnt, 12+i_pnt), 2, BRIDLE)

            Segment(9+i_seg, (9+i_pnt, 14+i_pnt), 2, BRIDLE)
            Segment(10+i_seg, (10+i_pnt, 13+i_pnt), 1, BRIDLE)
            Segment(11+i_seg, (11+i_pnt, 15+i_pnt), 2, BRIDLE)
            Segment(12+i_seg, (12+i_pnt, 17+i_pnt), 5, BRIDLE)
            
            Segment(13+i_seg, (13+i_pnt, 14+i_pnt), 1, BRIDLE)
            Segment(14+i_seg, (13+i_pnt, 15+i_pnt), 1, BRIDLE)
            
            Segment(15+i_seg, (14+i_pnt, 16+i_pnt), 3, BRIDLE)
            Segment(16+i_seg, (15+i_pnt, 16+i_pnt), 3, BRIDLE)
            Segment(17+i_seg, (15+i_pnt, 17+i_pnt), 3, BRIDLE)
        ]
        pulleys = [
            pulleys
            Pulley(1+i_pul, (13+i_seg, 14+i_seg), nothing)
            Pulley(2+i_pul, (16+i_seg, 17+i_seg), nothing)
        ]
        push!(attach_points, points[end-1])
        push!(attach_points, points[end])
        return nothing
    end

    function create_tether(attach_point, winch, type)
        l0 = s.set.l_tether / s.set.segments
        segment_idxs = Int16[]
        for i in 1:s.set.segments
            pos = attach_point.pos .+ [0, 0, -i*l0]
            i_pnt = length(points) # last point idx
            i_seg = length(segments) # last segment idx
            if i == s.set.segments
                points = [points; WinchPoint(1+i_pnt, pos, winch)]
                segments = [segments; Segment(1+i_seg, (i_pnt, 1+i_pnt), l0, type)]
            elseif i == 1
                points = [points; Point(1+i_pnt, pos)]
                segments = [segments; Segment(1+i_seg, (attach_point.idx, 1+i_pnt), l0, type)]
            else
                points = [points; Point(1+i_pnt, pos)]
                segments = [segments; Segment(1+i_seg, (i_pnt, 1+i_pnt), l0, type)]
            end
            push!(segment_idxs, 1+i_seg)
            i_pnt = length(points)
        end
        i_tether = length(tethers)
        winch_point = points[end].idx
        tethers = [tethers; Tether(1+i_tether, segment_idxs, winch_point)]
        return nothing
    end

    create_bridle(bridle_gammas[1:2], bridle_limits[1:2])
    create_bridle(bridle_gammas[3:4], bridle_limits[1:2])

    winches = [TorqueControlledMachine(s.set) for i in 1:4]
    create_tether.(attach_points, winches, [POWER, STEERING, POWER, STEERING])

    system = PointMassSystem(points, groups, segments, pulleys, tethers)
    s.point_system = system
    # plot(system, 0.0)
    return system
end

function create_sys!(s::KPSQ; init=false)
    system = create_point_mass_system!(s, s.wing)
    points, groups, segments, pulleys, tethers = 
        system.points, system.groups, system.segments, system.pulleys, system.tethers

    eqs = []
    tether_kite_force = zeros(3)
    tether_kite_torque = zeros(3)

    @parameters begin
        measured_wind_dir_gnd = s.measure.wind_dir_gnd
        measured_sphere_pos[1:2, 1:2] = s.measure.sphere_pos
        measured_sphere_vel[1:2, 1:2] = s.measure.sphere_vel
        measured_sphere_acc[1:2, 1:2] = s.measure.sphere_acc
        measured_tether_length[1:3] = s.measure.tether_length
        measured_tether_vel[1:3]    = s.measure.tether_vel
        measured_tether_acc[1:3]    = s.measure.tether_acc
    end
    @variables begin
        set_values(t)[1:3] # left right middle

        # potential differential variables
        kite_pos(t)[1:3] # xyz pos of kite in world frame
        kite_vel(t)[1:3]
        kite_acc(t)[1:3]
        distance(t)
        distance_vel(t)
        distance_acc(t)
        wind_scale_gnd(t)
        ω_p(t)[1:3] # turn rate in principal frame
        ω_b(t)[1:3] # turn rate in body frame
        α_p(t)[1:3] # angular acceleration in principal frame
        α_b(t)[1:3] # angular acceleration in body frame
        trailing_edge_angle(t)[eachindex(groups)] # angle left / right
        trailing_edge_ω(t)[eachindex(groups)] # angular rate
        trailing_edge_α(t)[eachindex(groups)] # angular acc
        twist_angle(t)[eachindex(groups)] # main body angle left / right
        twist_ω(t)[eachindex(groups)] # angular rate
        twist_α(t)[eachindex(groups)] # angular acc

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
        winch_force(t)[eachindex(tethers)] # normalized tether forces at the winch
        rho_kite(t)
        wind_vec_gnd(t)[1:3]
        wind_vel_kite(t)[1:3]
        va_kite(t)[1:3]
        va_kite_b(t)[1:3]
    end

    function force_eqs!()

        # ==================== POINTS ==================== #
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
                force[:, point.idx]  ~ F
            ]

            if point isa WinchPoint
                eqs = [
                    eqs
                    pos[:, point.idx]    ~ zeros(3)
                    vel[:, point.idx]    ~ zeros(3)
                    acc[:, point.idx]    ~ zeros(3)
                ]
            elseif point isa KitePoint
                tether_kite_force .+= F
                tether_kite_torque .+= (s.R_b_p * point.pos) × (R_p_w' * F)
                chord_b = point.pos - point.fixed_pos
                idx = point.pos[2] > 0 ? 1 : 2
                pos_b = point.fixed_pos + rotate_v_around_k(chord_b, point.y_panel, twist_angle[idx])
                pos_w = kite_pos + R_b_w * pos_b
                eqs = [
                    eqs
                    pos[:, point.idx]    ~ pos_w
                    vel[:, point.idx]    ~ zeros(3)
                    acc[:, point.idx]    ~ zeros(3)
                ]
            elseif point isa Point
                eqs = [
                    eqs
                    D(pos[:, point.idx]) ~ vel[:, point.idx]
                    D(vel[:, point.idx]) ~ acc[:, point.idx]
                    acc[:, point.idx]    ~ force[:, point.idx] / mass + G_EARTH
                ]
            else
                throw(ArgumentError("Unknown point type: $(typeof(point))"))
            end
        end

        # ==================== GROUPS ==================== #
        for group in groups
            eqs = [
                eqs
                D(twist_angle[group.idx]) ~ twist_ω[group.idx]
                D(twist_ω[group.idx]) ~ twist_α[group.idx]
                # twist_α[group.idx] ~ 
            ]
        end

        # ==================== SEGMENTS ==================== #
        @variables begin
            segment(t)[1:3, eachindex(segments)]
            unit_vector(t)[1:3, eachindex(segments)]
            len(t)[eachindex(segments)]
            l0(t)[eachindex(segments)]
            rel_vel(t)[1:3, eachindex(segments)]
            spring_vel(t)[eachindex(segments)]
            spring_force(t)[eachindex(segments)]
            spring_force_vec(t)[1:3, eachindex(segments)]
        end
        for segment in segments
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
                        eqs = [
                            eqs
                            l0[segment.idx] ~ tether_length[tether.idx] / length(tether.segments)
                        ]
                        in_tether += 1
                    end
                end
                (in_tether != 1) && throw(ArgumentError("Segment number $(segment.idx) is part of 
                    $in_tether tethers, and should be part of exactly 1 tether."))
            else
                throw(ArgumentError("Unknown segment type: $(segment.type)"))
            end

            p1, p2 = segment.points[1], segment.points[2]

            (segment.type === BRIDLE) && (diameter = s.bridle_tether_diameter)
            (segment.type === POWER) && (diameter = s.power_tether_diameter)
            (segment.type === STEERING) && (diameter = s.steering_tether_diameter)

            stiffness = s.set.e_tether * (diameter/2000)^2 * pi
            (segment.type === BRIDLE) && (compression_frac = 1.0)
            (segment.type === POWER) && (compression_frac = 0.1)
            (segment.type === STEERING) && (compression_frac = 0.1)
            
            damping = (s.set.damping / s.set.c_spring) * stiffness
            @show damping
            (segment.type === BRIDLE) && (damping = 10damping)

            eqs = [
                eqs
                # spring force equations
                segment[:, segment.idx]      ~ pos[:, p2] - pos[:, p1]
                len[segment.idx]             ~ norm(segment[:, segment.idx])
                unit_vector[:, segment.idx]  ~ segment[:, segment.idx]/len[segment.idx]
                rel_vel[:, segment.idx]      ~ vel[:, p1] - vel[:, p2]
                spring_vel[segment.idx]      ~ rel_vel[:, segment.idx] ⋅ unit_vector[:, segment.idx]
                spring_force[segment.idx]    ~ (stiffness * segment.l0 * (len[segment.idx] - l0[segment.idx]) - 
                                            damping * segment.l0 * spring_vel[segment.idx])
                stiffness[segment.idx]       ~ ifelse(len[segment.idx] > l0[segment.idx],
                                            stiffness / len[segment.idx],
                                            compression_frac * stiffness / len[segment.idx]
                )
                damping[segment.idx]         ~ damping / len[segment.idx]
                spring_force ~  (stiffness[segment.idx] * segment.l0 * (len[segment.idx] - l0[segment.idx]) - 
                                damping[segment.idx] * segment.l0 * spring_vel[segment.idx])
                spring_force_vec[:, segment.idx]  ~ spring_force[segment.idx] * unit_vector[:, segment.idx]

                # drag force equations
                height[segment.idx]          ~ max(0.0, 0.5(pos[:, p1][3] + pos[:, p2][3]))
                segment_vel[:, segment.idx]   ~ 0.5(vel[:, p1] + vel[:, p2])
                segment_rho[segment.idx]      ~ calc_rho(s.am, height[i])
                wind_vel[:, segment.idx]     ~ AtmosphericModels.calc_wind_factor(s.am, height[i], s.set.profile_law) * wind_vec_gnd
                va[:, segment.idx] ~ wind_vel - segment_vel
                area[segment.idx]            ~ len[segment.idx] * segment.diameter
                app_perp_vel[:, segment.idx] ~ va[:, segment.idx] - 
                                            (va[:, segment.idx] ⋅ unit_vector[:, segment.idx]) * unit_vector[:, segment.idx]
                drag_force                  ~ (0.5 * tether_rho[segment.idx] * s.set.cd_tether * norm(va[:, segment.idx]) * 
                                            area[segment.idx]) * app_perp_vel[:, segment.idx]
            ]
        end

        # ==================== PULLEYS ==================== #
        pulley_damping = 10
        for pulley in pulleys
            segment1, segment2 = segments[pulley.segments[1]], segments[pulley.segments[2]]
            pulley.sum_length = segment1.l0 + segment2.l0
            mass = pulley.sum_length * segment1.mass_per_meter
            eqs = [
                eqs
                D(pulley_l0) ~ pulley_vel
                D(pulley_vel) ~ pulley_acc
                pulley_force[pulley.idx]    ~ spring_force[pulley.segments[1]] - spring_force[pulley.segments[2]]
                pulley_acc[pulley.idx]      ~ pulley_force[pulley.idx] / mass - pulley_damping * pulley_vel[pulley.idx]
            ]
        end

        # ==================== TETHERS ==================== #
        @variables begin
            tether_length(t)[eachindex(tethers)]
            tether_vel(t)[eachindex(tethers)]
            tether_acc(t)[eachindex(tethers)]
            winch_force(t)[eachindex(tethers)]
        end
        for tether in tethers
            winch_point = nothing
            winch_found = 0
            for segment in tether.segments
                for idx in segments[segment].points
                    if points[idx] isa WinchPoint
                        winch_point = points[idx]
                        winch_found += 1
                    end
                end
            end
            (winch_found != 1) && throw(ArgumentError("Tether number $(tether.idx) has
                $winch_found winches, but should have exactly 1."))
            eqs = [
                eqs
                D(tether_length[tether.idx]) ~ tether_vel[tether.idx]
                D(tether_vel[tether.idx]) ~ tether_acc[tether.idx]
                tether_acc[tether.idx] ~ calc_torque_acc( # TODO: torque and speed control
                    winch_point.winch, tether_vel[tether.idx], 
                    winch_force[tether.idx], 
                    set_values[tether.idx]
                )
                winch_force[tether.idx] ~ norm(force[:, winch_point.idx])
            ]
        end
    end

    function diff_eqs!()
        Q_b_p = quaternion_conjugate(s.Q_p_b)
        if !init
            Ω = [0       -ω_p[1]  -ω_p[2]  -ω_p[3];
                ω_p[1]   0        ω_p[3]  -ω_p[2];
                ω_p[2]  -ω_p[3]   0        ω_p[1];
                ω_p[3]   ω_p[2]  -ω_p[1]   0]
            eqs = [
                eqs
                [D(Q_p_w[i]) ~ Q_vel[i] for i in 1:4]
                [Q_vel[i] ~ 0.5 * sum(Ω[i, j] * Q_p_w[j] for j in 1:4) for i in 1:4]
                [R_b_w[:, i] ~ quaternion_to_rotation_matrix(quaternion_multiply(Q_p_w, Q_b_p))[:, i] for i in 1:3] # https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation#Performance_comparisons
                [R_p_w[:, i] ~ quaternion_to_rotation_matrix(Q_p_w)[:, i] for i in 1:3]
                D(ω_p) ~ α_p
                ω_b ~ s.R_b_p' * ω_p
                α_p[1] ~ (torque_p[1] + (s.I_p[2] - s.I_p[3]) * ω_p[2] * ω_p[3]) / s.I_p[1]
                α_p[2] ~ (torque_p[2] + (s.I_p[3] - s.I_p[1]) * ω_p[3] * ω_p[1]) / s.I_p[2]
                α_p[3] ~ (torque_p[3] + (s.I_p[1] - s.I_p[2]) * ω_p[1] * ω_p[2]) / s.I_p[3]
                α_b ~ s.R_b_p' * α_p
                
                [D(kite_pos[i]) ~ kite_vel[i] for i in 1:3]
                [D(kite_vel[i]) ~ kite_acc[i] for i in 1:3]
                aero_kite_force ~ (R_b_w * s.vsm_solver.sol.force_coefficients) * (0.5 * rho_kite * norm(va_kite)^2) * s.aero.projected_area
                kite_acc        ~ (tether_kite_force + aero_kite_force) / s.set.mass

                distance            ~ norm(kite_pos)
                distance_vel        ~ kite_vel ⋅ normalize(kite_pos)
                distance_acc        ~ kite_acc ⋅ normalize(kite_pos)    
                D(wind_scale_gnd) ~ 0
            ]
        else
            idamp = 10
            scale = 1e-2

            # no movement around body z axis
            Ω = [0       -ω_b[1]  -ω_b[2]  -0     ;
                ω_b[1]   0        0       -ω_b[2];
                ω_b[2]  -0        0        ω_b[1];
                0        ω_b[2]  -ω_b[1]   0     ]
            
            # from measurements
            elevation     = mean(measured_sphere_pos[1, :])
            azimuth       = mean(measured_sphere_pos[2, :])
            elevation_vel = mean(measured_sphere_vel[1, :])
            azimuth_vel   = mean(measured_sphere_vel[2, :])
            elevation_acc = mean(measured_sphere_acc[1, :])
            azimuth_acc   = mean(measured_sphere_acc[2, :])

            r = (measured_sphere_pos[:, 2] - measured_sphere_pos[:, 1]) / 2
            perp_r = [-r[2], r[1]]
            rot_vel = (measured_sphere_vel[:, 1] - measured_sphere_vel[:, 2]) ⋅ (perp_r / norm(r))
            rot_acc = (measured_sphere_acc[:, 1] - measured_sphere_acc[:, 2]) ⋅ (perp_r / norm(r))
            ω_z = rot_vel / norm(r)
            α_z = rot_acc / norm(r)

            angular_acc = measured_tether_acc / s.set.drum_radius
            net_torque = angular_acc * s.set.inertia_total
            measured_winch_force = (net_torque - set_values) / s.set.drum_radius

            # scaled down, partially fixed, and damped equations to find unmeasured variables
            eqs = [
                eqs
                [D(Q_b_w[i]) ~ Q_vel[i] for i in 1:4]
                Q_p_w ~ quaternion_multiply(Q_b_w, quaternion_conjugate(Q_b_p))
                [Q_vel[i] ~ 0.5 * sum(Ω[i, j] * Q_b_w[j] for j in 1:4) for i in 1:4]
                [R_b_w[:, i] ~ quaternion_to_rotation_matrix(Q_b_w)[:, i] for i in 1:3] # https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation#Performance_comparisons
                [R_p_w[:, i] ~ quaternion_to_rotation_matrix(Q_p_w)[:, i] for i in 1:3]
                D(ω_b[1:2]) ~ scale * α_b[1:2] - idamp * ω_b[1:2]
                ω_b[3] ~ ω_z
                ω_p ~ s.R_b_p * ω_b
                torque_b ~ s.R_b_p' * torque_p
                α_b[1] ~ (torque_b[1]) / s.I_b[1]
                α_b[2] ~ (torque_b[2]) / s.I_b[2]
                α_b[3] ~ α_z

                kite_pos        ~ rotate_around_z(rotate_around_y([distance, 0, 0], -elevation), azimuth)
                kite_vel        ~ distance_vel * normalize(kite_pos - pos[:, 3]) + 
                                    rotate_around_z(rotate_around_y([0, azimuth_vel * distance, elevation_vel * distance], -elevation), azimuth)
                kite_acc        ~ distance_acc * normalize(kite_pos - pos[:, 3]) + 
                                    rotate_around_z(rotate_around_y([0, azimuth_acc * distance, elevation_acc * distance], -elevation), azimuth)
                D(distance)     ~ distance_vel
                D(distance_vel) ~ distance_acc - idamp * distance_vel
                distance_acc    ~ scale * (measured_tether_acc[3] - tether_acc[3])
                D(wind_scale_gnd)  ~ 100scale * (mean(measured_winch_force[1:2]) - mean(winch_force[1:2]))

                [pos[:, i]              .~ 0.0 for i in 1:3]
                [D.(pos[:, i])          .~ vel[:, i] for i in 4:s.i_A-1]
                D(trailing_edge_angle)   ~ trailing_edge_ω
                [vel[:, i]              .~ 0.0 for i in 1:3]
                [D.(vel[:, i])          .~ scale * acc[:, i] - idamp * vel[:, i] for i in 4:s.i_A-1] # TODO: ADD CENTRIFUGAL FORCE DUE TO KITE ROTATION
                D(trailing_edge_ω)       ~ 0.1scale * trailing_edge_α - idamp * trailing_edge_ω
                tether_length           ~ measured_tether_length
                tether_vel              ~ measured_tether_vel
            ]
        end
        return nothing
    end
    function scalar_eqs!()
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
            power_angle(t) # average flap angle
            power_vel(t)
            steering_angle(t) # difference between left and right flap angle
            steering_vel(t)
            tether_diff(t)
            tether_diff_vel(t)
            set_diff(t)
            azimuth(t)
            azimuth_vel(t)
            azimuth_acc(t)
            elevation(t)
            elevation_vel(t)
            elevation_acc(t)
            x_acc(t)
            y_acc(t)
            left_diff(t)
            right_diff(t)
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
            power_angle         ~ (trailing_edge_angle[1] + trailing_edge_angle[2]) / 2
            power_vel           ~ (trailing_edge_ω[1] + trailing_edge_ω[2]) / 2
            steering_angle      ~ trailing_edge_angle[2] - trailing_edge_angle[1]
            steering_vel        ~ trailing_edge_ω[2] - trailing_edge_ω[1]
            tether_diff         ~ tether_length[2] - tether_length[1]
            tether_diff_vel     ~ tether_vel[2] - tether_vel[1]
            set_diff            ~ set_values[2] - set_values[1]

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
            left_diff           ~ tether_length[1] - tether_length[3]
            right_diff          ~ tether_length[2] - tether_length[3]
        ]
        return nothing
    end

    force_eqs!()
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

    if !init
        discrete_events = [true => [Q_p_w[i] ~ normalize(Q_p_w)[i] for i in 1:4]]
    else
        discrete_events = [true => [Q_b_w[i] ~ normalize(Q_b_w)[i] for i in 1:4]]
    end
    
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
function model!(s::KPSQ; init=false)
    # init_pos!(s)
    
    sys, inputs = create_sys!(s; init)
    # structural_simplify(sys, (inputs, []))
    (sys, _) = structural_simplify(sys, (inputs, []); fully_determined=true)

    if !init
        u0map = [
            [sys.Q_p_w[i] => s.Q_p_w[i] for i in 1:4]
            [sys.ω_p[i] => 0 for i in 1:3]

            [sys.kite_pos[j] => s.kite_pos[j] for j in 1:3]
            [sys.kite_vel[i] => 0.0 for i in 1:3]

            [sys.pos[j, i] => s.pos[j, i] for j in 1:3 for i in 4:s.i_A-1]
            [sys.vel[j, i] => 0 for j in 1:3 for i in 4:s.i_A-1]

            [sys.trailing_edge_angle[i] => s.te_angle[i] for i in 1:2]
            [sys.trailing_edge_ω[i] => 0 for i in 1:2]

            [sys.tether_length[i] => s.measure.tether_length[i] for i in 1:3]
            [sys.tether_vel[j] => 0 for j in 1:3]

            sys.wind_scale_gnd => s.set.v_wind
        ]
    else
        Q_b_w = quaternion_multiply(s.Q_p_w, quaternion_conjugate(s.Q_p_b))
        u0map = [
            [sys.Q_b_w[i] => Q_b_w[i] for i in 1:4]
            [sys.ω_b[i] => 0 for i in 1:2]

            sys.distance => norm(s.kite_pos)
            sys.distance_vel => 0.0

            [sys.pos[j, i] => s.pos[j, i] for j in 1:3 for i in 4:s.i_A-1]
            [sys.vel[j, i] => 0 for j in 1:3 for i in 4:s.i_A-1]

            [sys.trailing_edge_angle[i] => s.te_angle[i] for i in 1:2]
            [sys.trailing_edge_ω[i] => 0 for i in 1:2]

            sys.wind_scale_gnd => s.set.v_wind
        ]
    end
    p0map = [sys.set_values[j] => s.measure.set_values[j] for j in 1:3]
    return sys, u0map, p0map
end
