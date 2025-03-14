
function init!(system::PointMassSystem, s::KPSQ)
    points, groups, segments, pulleys, tethers, winches = 
        system.points, system.groups, system.segments, system.pulleys, system.tethers, system.winches

    for segment in segments
        (segment.type === BRIDLE) && (segment.diameter = s.bridle_tether_diameter)
        (segment.type === POWER) && (segment.diameter = s.power_tether_diameter)
        (segment.type === STEERING) && (segment.diameter = s.steering_tether_diameter)
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
                tether_length += segment.l0
            end
        end
        winch.tether_length = tether_length
        @assert !(winch.tether_length ≈ 0)
    end

    min_pos = zeros(SimFloat, 3)
    min_pos[3] = Inf
    for point in points
        if point.pos[3] < min_pos[3]
            min_pos .= point.pos
        end
    end
    for point in points
        point.pos .-= min_pos
    end
    s.kite_pos .= -min_pos
    return nothing
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
    @show bridle_fracs
    for (i, gamma) in enumerate(bridle_gamma)
        for (j, frac) in enumerate(bridle_fracs)
            bridle_pos_b[:, j, i] .= s.pos_circle_center_b + 
                [s.leading_edge(gamma) + frac * (s.trailing_edge(gamma) - s.leading_edge(gamma)), cos(gamma) * s.set.radius, sin(gamma) * s.set.radius]
        end
    end
    return bridle_pos_b
end

function calc_inertia(s::KPSQ, wing::KiteWing)
    s.I_b .= [wing.inertia_tensor[1,1], wing.inertia_tensor[2,2], wing.inertia_tensor[3,3]]
    
    # Find principal axes
    eigenvals, eigenvecs = eigen(I)
    
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
    Q_p_b = rotation_matrix_to_quaternion(s.R_b_p')
    return I_p, R_b_p, Q_p_b
end

function measure_to_q(measure::Measurement, R_b_p)
    x = [1, 0, 0]
    z = [0, 0, 1]
    x = rotate_around_y(x, -measure.elevation)
    z = rotate_around_y(z, -measure.elevation)
    x = rotate_around_z(x, measure.azimuth)
    z = rotate_around_z(z, measure.azimuth)
    R_b_w = hcat(x, z × x, z)
    R_p_w = R_b_w * R_b_p'
    Q_p_w = rotation_matrix_to_quaternion(R_p_w)
    return Q_p_w
end

