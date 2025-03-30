
function PointMassSystem(s::KPSQ, wing::RamAirWing)
    # TODO: move as much of the code as possible from create_point_mass_system to other places, to make model creation easier.
    # 1. move bridle gamma calculation
    # 2. ...

    points = Point[]
    groups = KitePointGroup[]
    segments = Segment[]
    pulleys = Pulley[]
    tethers = Tether[]
    winches = Winch[]

    attach_points = Point[]
    
    bridle_top_left = [s.wing.R_cad_body * (s.top_bridle_points[i] + s.wing.T_cad_body) for i in eachindex(s.top_bridle_points)] # cad to kite frame
    bridle_top_right = [bridle_top_left[i] .* [1, -1, 1] for i in eachindex(s.top_bridle_points)]

    function create_bridle(bridle_top, gammas)
        i_pnt = length(points) # last point idx
        i_seg = length(segments) # last segment idx
        i_pul = length(pulleys) # last pulley idx

        i = 1
        for gamma in gammas # 2 gammas per bridle system
            le_pos = [wing.le_interp[i](gamma) for i in 1:3]
            chord = [wing.te_interp[i](gamma) for i in 1:3] .- le_pos
            y_airf = normalize([wing.le_interp[i](gamma-0.01) for i in 1:3] - le_pos)

            point_idxs = Int16[]
            for frac in s.bridle_fracs # 4 fracs
                pos = le_pos .+ chord .* frac
                points = [points; Point(i+i_pnt, pos, KITE)]
                push!(point_idxs, points[end].idx)
                i += 1
            end
            
            i_grp = 1 + length(groups)
            groups = [groups; KitePointGroup(i_grp, point_idxs, 1, chord, y_airf, DYNAMIC)]
        end

        points = [
            points
            Point(9+i_pnt, bridle_top[1], DYNAMIC)
            Point(10+i_pnt, bridle_top[2], DYNAMIC)
            Point(11+i_pnt, bridle_top[3], DYNAMIC)
            Point(12+i_pnt, bridle_top[4], DYNAMIC)

            Point(13+i_pnt, bridle_top[2] .+ [0, 0, -1], DYNAMIC)

            Point(14+i_pnt, bridle_top[1] .+ [0, 0, -2], DYNAMIC)
            Point(15+i_pnt, bridle_top[3] .+ [0, 0, -2], DYNAMIC)

            Point(16+i_pnt, bridle_top[1] .+ [0, 0, -3], DYNAMIC)
            Point(17+i_pnt, bridle_top[3] .+ [0, 0, -3], DYNAMIC)
        ]
        l1 = norm(points[9+i_pnt].pos_b - points[1+i_pnt].pos_b)
        l2 = norm(points[9+i_pnt].pos_b - points[5+i_pnt].pos_b)
        segments = [
            segments
            Segment(1+i_seg, (1+i_pnt, 9+i_pnt), BRIDLE, l1)
            Segment(2+i_seg, (2+i_pnt, 10+i_pnt), BRIDLE, l1)
            Segment(3+i_seg, (3+i_pnt, 11+i_pnt), BRIDLE, l1)
            Segment(4+i_seg, (4+i_pnt, 12+i_pnt), BRIDLE, l1)

            Segment(5+i_seg, (5+i_pnt, 9+i_pnt), BRIDLE, l2)
            Segment(6+i_seg, (6+i_pnt, 10+i_pnt), BRIDLE, l2)
            Segment(7+i_seg, (7+i_pnt, 11+i_pnt), BRIDLE, l2)
            Segment(8+i_seg, (8+i_pnt, 12+i_pnt), BRIDLE, l2)

            Segment(9+i_seg, (9+i_pnt, 14+i_pnt), BRIDLE, 2)
            Segment(10+i_seg, (10+i_pnt, 13+i_pnt), BRIDLE, 1)
            Segment(11+i_seg, (11+i_pnt, 15+i_pnt), BRIDLE, 2)
            Segment(12+i_seg, (12+i_pnt, 17+i_pnt), BRIDLE, 3)
            
            Segment(13+i_seg, (13+i_pnt, 14+i_pnt), BRIDLE, 1)
            Segment(14+i_seg, (13+i_pnt, 15+i_pnt), BRIDLE, 1)
            
            Segment(15+i_seg, (14+i_pnt, 16+i_pnt), BRIDLE, 1)
            Segment(16+i_seg, (15+i_pnt, 16+i_pnt), BRIDLE, 1)
            Segment(17+i_seg, (15+i_pnt, 17+i_pnt), BRIDLE, 1)
        ]
        pulleys = [
            pulleys
            Pulley(1+i_pul, (13+i_seg, 14+i_seg), DYNAMIC)
            Pulley(2+i_pul, (16+i_seg, 17+i_seg), DYNAMIC)
        ]
        push!(attach_points, points[end-1])
        push!(attach_points, points[end])
        return nothing
    end

    function create_tether(attach_point, type)
        l0 = s.set.l_tether / s.set.segments
        segment_idxs = Int16[]
        for i in 1:s.set.segments
            frac = i / s.set.segments
            pos = [(1-frac) * attach_point.pos_b[1], 
                    (1-frac) * attach_point.pos_b[2],
                    attach_point.pos_b[3] - i*l0]
            i_pnt = length(points) # last point idx
            i_seg = length(segments) # last segment idx
            if i == 1
                points = [points; Point(1+i_pnt, pos, DYNAMIC)]
                segments = [segments; Segment(1+i_seg, (attach_point.idx, 1+i_pnt), type)]
            elseif i == s.set.segments
                points = [points; Point(1+i_pnt, pos, WINCH)]
                segments = [segments; Segment(1+i_seg, (i_pnt, 1+i_pnt), type)]
            else
                points = [points; Point(1+i_pnt, pos, DYNAMIC)]
                segments = [segments; Segment(1+i_seg, (i_pnt, 1+i_pnt), type)]
            end
            push!(segment_idxs, 1+i_seg)
            i_pnt = length(points)
        end
        i_tether = length(tethers)
        winch_point_idx = points[end].idx
        tethers = [tethers; Tether(1+i_tether, segment_idxs, winch_point_idx)]
        return tethers[end].idx
    end

    gammas = [-3/4, -1/4, 1/4, 3/4] * wing.gamma_tip
    create_bridle(bridle_top_left, gammas[1:2])
    create_bridle(bridle_top_right, gammas[3:4])

    left_power_idx = create_tether(attach_points[1], POWER)
    right_power_idx = create_tether(attach_points[3], POWER)
    left_steering_idx = create_tether(attach_points[2], STEERING)
    right_steering_idx = create_tether(attach_points[4], STEERING)

    winches = [winches; Winch(1, TorqueControlledMachine(s.set), [left_power_idx, right_power_idx])]
    winches = [winches; Winch(2, TorqueControlledMachine(s.set), [left_steering_idx])]
    winches = [winches; Winch(3, TorqueControlledMachine(s.set), [right_steering_idx])]

    return PointMassSystem(points, groups, segments, pulleys, tethers, winches)
end


function init!(system::PointMassSystem, s::KPSQ, R_b_w)
    points, groups, segments, pulleys, tethers, winches = 
        system.points, system.groups, system.segments, system.pulleys, system.tethers, system.winches

    for segment in segments
        (segment.type === BRIDLE) && (segment.diameter = 0.001s.bridle_tether_diameter)
        (segment.type === POWER) && (segment.diameter = 0.001s.power_tether_diameter)
        (segment.type === STEERING) && (segment.diameter = 0.001s.steering_tether_diameter)
        (segment.l0 ≈ 0) && (segment.l0 = norm(points[segment.points[1]].pos_b - points[segment.points[2]].pos_b) * 0.9999)
        @assert (0 < segment.diameter < 1)
        @assert (segment.l0 > 0)
    end

    for pulley in pulleys
        segment1, segment2 = segments[pulley.segments[1]], segments[pulley.segments[2]]
        pulley.sum_length = segment1.l0 + segment2.l0
        @assert !(pulley.sum_length ≈ 0)
    end

    for winch in winches
        tether_length = 0.0
        for tether in tethers[winch.tethers]
            for segment in segments[tether.segments]
                tether_length += segment.l0 / length(winch.tethers) * 0.9999
            end
        end
        winch.tether_length = tether_length
        @assert !(winch.tether_length ≈ 0)
    end

    min_z = Inf
    for point in points
        if point.pos_b[3] < min_z
            min_z = point.pos_b[3]
        end
    end
    for point in points
        point.pos_w[3] -= min_z
        point.pos_w .= R_b_w * point.pos_w
    end
    init_kite_pos = R_b_w * [0.0, 0.0, -min_z]
    return init_kite_pos
end

function measure_to_q(measure::Measurement)
    x = [0, 0, -1] # laying flat along x axis
    z = [1, 0, 0] # laying flat along x axis
    x = rotate_around_y(x, -measure.elevation)
    z = rotate_around_y(z, -measure.elevation)
    x = rotate_around_z(x, measure.azimuth)
    z = rotate_around_z(z, measure.azimuth)
    R_b_w = hcat(x, z × x, z)
    Q_b_w = rotation_matrix_to_quaternion(R_b_w)
    return Q_b_w, R_b_w
end

