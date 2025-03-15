
function PointMassSystem(s::KPSQ, wing::KiteWing)
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
    
    bridle_gammas, bridle_limits = find_bridle_gammas!(s, wing)

    function create_bridle(gammas, limits)
        i_pnt = length(points) # last point idx
        i_seg = length(segments) # last segment idx
        i_pul = length(pulleys) # last pulley idx

        i = 1
        for (gamma, limit) in zip(gammas, limits) # 2 gammas with 2 pairs of limits
            le_pos = [wing.le_interp[i](gamma) for i in 1:3]
            chord = [wing.te_interp[i](gamma) for i in 1:3] .- le_pos
            fixed_pos = le_pos .+ chord .* s.bridle_fracs[1]
            y_panel = normalize([wing.le_interp[i](gamma-0.01) for i in 1:3] - le_pos)

            point_idxs = Int16[]
            for frac in s.bridle_fracs # 4 fracs
                pos = le_pos .+ chord .* frac
                points = [points; Point(i+i_pnt, pos, KITE)]
                push!(point_idxs, points[end].idx)
                i += 1
            end
            
            i_grp = 1 + length(groups)
            y_lim = (wing.le_interp[2](limit[1]), wing.le_interp[2](limit[2])) # TODO: ylim is slightly off-centre
            groups = [groups; KitePointGroup(i_grp, point_idxs, y_lim, fixed_pos, chord, y_panel)]
        end

        mean_le = [wing.le_interp[i](mean(gammas)) for i in 1:3]
        chord_length = norm([wing.te_interp[i](mean(gammas)) for i in 1:3] .- mean_le)
        xs = s.bridle_fracs .* chord_length
        bridle_top = mean_le .+ [0, 0, -3]

        points = [
            points
            Point(9+i_pnt, bridle_top .+ [xs[1], 0, 0], DYNAMIC)
            Point(10+i_pnt, bridle_top .+ [xs[2], 0, 0], DYNAMIC)
            Point(11+i_pnt, bridle_top .+ [xs[3], 0, 0], DYNAMIC)
            Point(12+i_pnt, bridle_top .+ [xs[4], 0, 0], DYNAMIC)

            Point(13+i_pnt, bridle_top .+ [xs[2], 0, -1], DYNAMIC)

            Point(14+i_pnt, bridle_top .+ [xs[1], 0, -2], DYNAMIC)
            Point(15+i_pnt, bridle_top .+ [xs[3], 0, -2], DYNAMIC)

            Point(16+i_pnt, bridle_top .+ [xs[1], 0, -5], DYNAMIC)
            Point(17+i_pnt, bridle_top .+ [xs[3], 0, -5], DYNAMIC)
        ]
        segments = [
            segments
            Segment(1+i_seg, (1+i_pnt, 9+i_pnt), BRIDLE)
            Segment(2+i_seg, (2+i_pnt, 10+i_pnt), BRIDLE)
            Segment(3+i_seg, (3+i_pnt, 11+i_pnt), BRIDLE)
            Segment(4+i_seg, (4+i_pnt, 12+i_pnt), BRIDLE)

            Segment(5+i_seg, (5+i_pnt, 9+i_pnt), BRIDLE)
            Segment(6+i_seg, (6+i_pnt, 10+i_pnt), BRIDLE)
            Segment(7+i_seg, (7+i_pnt, 11+i_pnt), BRIDLE)
            Segment(8+i_seg, (8+i_pnt, 12+i_pnt), BRIDLE)

            Segment(9+i_seg, (9+i_pnt, 14+i_pnt), BRIDLE)
            Segment(10+i_seg, (10+i_pnt, 13+i_pnt), BRIDLE)
            Segment(11+i_seg, (11+i_pnt, 15+i_pnt), BRIDLE)
            Segment(12+i_seg, (12+i_pnt, 17+i_pnt), BRIDLE)
            
            Segment(13+i_seg, (13+i_pnt, 14+i_pnt), BRIDLE)
            Segment(14+i_seg, (13+i_pnt, 15+i_pnt), BRIDLE)
            
            Segment(15+i_seg, (14+i_pnt, 16+i_pnt), BRIDLE)
            Segment(16+i_seg, (15+i_pnt, 16+i_pnt), BRIDLE)
            Segment(17+i_seg, (15+i_pnt, 17+i_pnt), BRIDLE)
        ]
        pulleys = [
            pulleys
            Pulley(1+i_pul, (13+i_seg, 14+i_seg))
            Pulley(2+i_pul, (16+i_seg, 17+i_seg))
        ]
        push!(attach_points, points[end-1])
        push!(attach_points, points[end])
        return nothing
    end

    function create_tether(attach_point, type)
        l0 = s.set.l_tether / s.set.segments
        segment_idxs = Int16[]
        for i in 1:s.set.segments
            pos = attach_point.pos_b .+ [0, 0, -i*l0]
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

    create_bridle(bridle_gammas[1:2], bridle_limits[1:2])
    create_bridle(bridle_gammas[3:4], bridle_limits[3:4])

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
        (segment.type === BRIDLE) && (segment.diameter = s.bridle_tether_diameter)
        (segment.type === POWER) && (segment.diameter = s.power_tether_diameter)
        (segment.type === STEERING) && (segment.diameter = s.steering_tether_diameter)
        segment.l0 = norm(points[segment.points[1]].pos_b - points[segment.points[2]].pos_b)
        @assert !(segment.diameter ≈ 0)
        @assert !(segment.l0 ≈ 0)
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
                tether_length += segment.l0 / length(winch.tethers)
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
        point.pos_w .= R_b_w * point.pos_b
    end
    init_kite_pos = R_b_w * [0.0, 0.0, -min_z]
    return init_kite_pos
end


function find_bridle_gammas!(s::KPSQ, wing::KiteWing; n_groups=4)
    @assert iseven(n_groups) "Number of groups must be even"
    half_area = wing.area_interp(0.0)
    
    # Create target areas - evenly distributed on left side of the wing
    target_areas = [(i / n_groups) * (half_area) for i in 1:n_groups-1]
    
    function equations!(F, gamma, p)
        for i in 1:n_groups-1
            F[i] = wing.area_interp(gamma[i]) - target_areas[i]
        end
    end
    
    # Initial guess: evenly spaced between -gamma_tip and middle
    gamma_tip = wing.gamma_tip
    gamma0 = [-gamma_tip + (i / n_groups) * gamma_tip for i in 1:n_groups-1]
    
    prob = NonlinearProblem(equations!, gamma0, nothing)
    result = NonlinearSolve.solve(prob, NewtonRaphson())
    
    # Mirror the solution for both sides
    bridle_gamma = zeros(n_groups)
    bridle_gamma[1:n_groups÷2] .= result.u[1:2:end]
    bridle_gamma[n_groups÷2+1:end] .= -reverse(result.u[1:2:end])

    limits = [zeros(2) for _ in 1:n_groups]
    for i in eachindex(limits[1:n_groups÷2])
        if i == 1
            limits[i] .= (-gamma_tip, result.u[2i])
        elseif 2i == length(result.u)+1
            limits[i] .= (result.u[2i-2], 0.0)
        else
            limits[i] .= (result.u[2i-2], result.u[2i])
        end
        limits[n_groups+1-i] .= -reverse(limits[i])
    end
    return bridle_gamma, limits
end


function init_bridle_pos!(s::KPSQ)
    bridle_pos_b = zeros(SimFloat, 3, 4, 4) # xyz, length, width
    bridle_gamma = zeros(SimFloat, 4)

    calc_inertia!(s, s.wing)
    find_attachment_gamma!(s, bridle_gamma)

    bridle_fracs = (s.set.bridle_connect[2:5] .- s.set.bridle_connect[1]) / (s.set.bridle_connect[1] - s.set.bridle_connect[5])
    for (i, gamma) in enumerate(bridle_gamma)
        for (j, frac) in enumerate(bridle_fracs)
            bridle_pos_b[:, j, i] .= s.pos_circle_center_b + 
                [s.leading_edge(gamma) + frac * (s.trailing_edge(gamma) - s.leading_edge(gamma)), cos(gamma) * s.set.radius, sin(gamma) * s.set.radius]
        end
    end
    return bridle_pos_b
end

function calc_inertia(wing::KiteWing)
    I_b = [wing.inertia_tensor[1,1], wing.inertia_tensor[2,2], wing.inertia_tensor[3,3]]
    
    # Find principal axes
    eigenvals, eigenvecs = eigen(wing.inertia_tensor)
    
    # Sort by magnitude
    p = sortperm(eigenvals)
    eigenvals = eigenvals[p]
    eigenvecs = eigenvecs[:, p]
    
    # Ensure right-handed coordinate system
    if det(eigenvecs) < 0
        eigenvecs[:, 3] .*= -1
    end

    # Store results
    I_p = eigenvals
    R_b_p = eigenvecs
    Q_p_b = rotation_matrix_to_quaternion(R_b_p')
    return I_p, I_b, R_b_p, Q_p_b
end

function measure_to_q(measure::Measurement, R_b_p)
    x = [0, 0, -1] # laying flat along x axis
    z = [1, 0, 0] # laying flat along x axis
    x = rotate_around_y(x, -measure.elevation)
    z = rotate_around_y(z, -measure.elevation)
    x = rotate_around_z(x, measure.azimuth)
    z = rotate_around_z(z, measure.azimuth)
    R_b_w = hcat(x, z × x, z)
    R_p_w = R_b_w * R_b_p'
    Q_p_w = rotation_matrix_to_quaternion(R_p_w)
    return Q_p_w, R_b_w
end

