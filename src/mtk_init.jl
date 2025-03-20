
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
    
    bridle_gammas, bridle_limits = find_bridle_gammas!(s, wing)
    bridle_top_left = [s.top_bridle_points[i] - s.wing.center_of_mass for i in eachindex(s.top_bridle_points)] # cad to kite frame
    bridle_top_right = [bridle_top_left[i] .* [1, -1, 1] for i in eachindex(s.top_bridle_points)]

    function create_bridle(gammas, limits, bridle_top)
        i_pnt = length(points) # last point idx
        i_seg = length(segments) # last segment idx
        i_pul = length(pulleys) # last pulley idx

        i = 1
        for (gamma, limit) in zip(gammas, limits) # 2 gammas with 2 pairs of limits
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
            y_lim = (wing.le_interp[2](limit[1]), wing.le_interp[2](limit[2])) # TODO: ylim is slightly off-centre
            groups = [groups; KitePointGroup(i_grp, point_idxs, y_lim, 1, chord, y_airf, DYNAMIC)]
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

    create_bridle(bridle_gammas[1:2], bridle_limits[1:2], bridle_top_left)
    create_bridle(bridle_gammas[3:4], bridle_limits[3:4], bridle_top_right)

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


function find_bridle_gammas!(s::KPSQ, wing::RamAirWing; n_groups=4)
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

function calc_inertia(I_b_tensor)
    I_b = [I_b_tensor[1,1], I_b_tensor[2,2], I_b_tensor[3,3]]
    
    # Find principal axes
    eigenvals, eigenvecs = eigen(I_b_tensor)
    
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

function calc_inertia2(I_b_tensor)
    
    # Find principal axes (eigenvalues and eigenvectors)
    eigenvals, eigenvecs = eigen(I_b_tensor)
    
    # Sort by magnitude of eigenvalues
    p = sortperm(eigenvals)
    I_p = eigenvals[p]
    eigenvecs = eigenvecs[:, p]
    
    # Ensure right-handed coordinate system
    if det(eigenvecs) < 0
        eigenvecs[:, 3] *= -1
    end
    
    # Find the best alignment with the standard basis [1,0,0], [0,1,0], [0,0,1]
    # Calculate how well each eigenvector aligns with each standard basis vector
    alignment = abs.([
        eigenvecs[:,1]⋅[1,0,0] eigenvecs[:,1]⋅[0,1,0] eigenvecs[:,1]⋅[0,0,1];
        eigenvecs[:,2]⋅[1,0,0] eigenvecs[:,2]⋅[0,1,0] eigenvecs[:,2]⋅[0,0,1];
        eigenvecs[:,3]⋅[1,0,0] eigenvecs[:,3]⋅[0,1,0] eigenvecs[:,3]⋅[0,0,1]
    ])
    
    # Find best assignment of eigenvectors to coordinate axes
    assigned_cols = zeros(Int, 3)
    assigned_rows = zeros(Int, 3)
    
    # Greedy algorithm to find best matching
    for _ in 1:3
        # Find max value in alignment matrix among unassigned rows/columns
        max_val = -1.0
        best_row = 0
        best_col = 0
        
        for row in 1:3
            if assigned_rows[row] == 0
                for col in 1:3
                    if assigned_cols[col] == 0 && alignment[row, col] > max_val
                        max_val = alignment[row, col]
                        best_row = row
                        best_col = col
                    end
                end
            end
        end
        
        # Assign this eigenvector to this axis
        assigned_rows[best_row] = best_col
        assigned_cols[best_col] = best_row
    end
    
    # Create reordered eigenvector matrix and reordered eigenvalues
    R_b_p = zeros(3, 3)
    I_p_reordered = zeros(3)
    
    for i in 1:3
        axis = assigned_rows[i]
        # Make sure the eigenvector points in the positive direction of standard basis
        sign_factor = sign(eigenvecs[:,i]⋅[axis==1, axis==2, axis==3])
        if sign_factor == 0
            sign_factor = 1  # Default to positive if dot product is zero
        end
        
        R_b_p[:,axis] = sign_factor * eigenvecs[:,i]
        I_p_reordered[axis] = I_p[i]
    end
    
    # Ensure the matrix is still a rotation matrix (det = 1)
    if det(R_b_p) < 0
        R_b_p[:,3] *= -1
    end
    
    # Calculate the quaternion representation
    Q_p_b = rotation_matrix_to_quaternion(R_b_p')
    
    # Get the diagonal terms of the original tensor 
    I_b = [I_b_tensor[1,1], I_b_tensor[2,2], I_b_tensor[3,3]]
    
    return I_p_reordered, I_b, R_b_p, Q_p_b
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

