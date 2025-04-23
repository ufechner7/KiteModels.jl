
function VortexStepMethod.RamAirWing(set::Settings; prn=true, kwargs...)
    obj_path = joinpath(dirname(get_data_path()), set.model)
    dat_path = joinpath(dirname(get_data_path()), set.foil_file)
    return RamAirWing(obj_path, dat_path; 
        mass=set.mass, crease_frac=set.crease_frac, n_groups=length(set.bridle_fracs), 
        align_to_principal=true, prn, kwargs...
    )
end


@enum SegmentType begin
    POWER
    STEERING
    BRIDLE
end

@enum DynamicsType begin
    DYNAMIC
    STATIC
    KITE
    WINCH
end

"""
A normal freely moving tether point
"""
mutable struct Point
    idx::Int16
    pos_b::KVec3 # pos relative to kite COM in body frame
    pos_w::KVec3 # pos in world frame
    type::DynamicsType
end
function Point(idx, pos_b, type)
    Point(idx, pos_b, copy(pos_b), type)
end

"""
Set of bridle lines that share the same twist angle and trailing edge angle and that rotates around the leading edge
"""
struct KitePointGroup
    idx::Int16
    points::Vector{Int16}
    le_point::KVec3 # point which the group rotates around under kite deformation
    chord::KVec3 # chord vector in body frame which the group rotates around under kite deformation
    y_airf::KVec3 # spanwise vector in local panel frame which the group rotates around under kite deformation
    type::DynamicsType
end

"""
A segment from one point index to another point index
"""
mutable struct Segment
    idx::Int16
    points::Tuple{Int16, Int16}
    type::SegmentType
    l0::SimFloat
    diameter::SimFloat
end
function Segment(idx, points, type)
    Segment(idx, points, type, zero(SimFloat), zero(SimFloat))
end
function Segment(idx, points, type, l0)
    Segment(idx, points, type, l0, zero(SimFloat))
end

"""
A pulley described by two segments with the common point of the segments being the pulley
"""
mutable struct Pulley
    idx::Int16
    segments::Tuple{Int16, Int16}
    type::DynamicsType
    sum_length::SimFloat
    function Pulley(idx, segments, type)
        new(idx, segments, type, zero(SimFloat))
    end
end

"""
A set of segments making a flexible tether. The winch point should only be part of one segment.
"""
struct Tether
    idx::Int16
    segments::Vector{Int16}
    winch_point::Int16
end

"""
A set of tethers or just one tether connected to a winch
"""
mutable struct Winch
    idx::Int16
    model::AbstractWinchModel
    tethers::Vector{Int16}
    tether_length::Float64
    function Winch(idx, model, tethers)
        new(idx, model, tethers, zero(Float64))
    end
end

struct PointMassSystem
    name::String
    points::Vector{Point}
    groups::Vector{KitePointGroup}
    segments::Vector{Segment}
    pulleys::Vector{Pulley}
    tethers::Vector{Tether}
    winches::Vector{Winch}
    function PointMassSystem(name, points, groups, segments, pulleys, tethers, winches)
        for (i, point) in enumerate(points)
            @assert point.idx == i
        end
        for (i, group) in enumerate(groups)
            @assert group.idx == i
        end
        for (i, segment) in enumerate(segments)
            @assert segment.idx == i
        end
        for (i, pulley) in enumerate(pulleys)
            @assert pulley.idx == i
        end
        for (i, tether) in enumerate(tethers)
            @assert tether.idx == i
        end
        for (i, winch) in enumerate(winches)
            @assert winch.idx == i
        end
        return new(name, points, groups, segments, pulleys, tethers, winches)
    end
end


function PointMassSystem(set::Settings, wing::RamAirWing)
    set.physical_model == "simple_ram" && length(set.bridle_fracs) != 2 && throw(ArgumentError("Model simple_ram should have 2 bridle fracs"))
    set.physical_model == "ram" && length(set.bridle_fracs) != 4 && throw(ArgumentError("Model ram should have 4 bridle fracs"))

    if set.physical_model == "ram"
        return create_ram_point_system(set, wing)
    elseif set.physical_model == "simple_ram"
        return create_simple_ram_point_system(set, wing)
    else
        throw(ArgumentError("Undefined physical model"))
    end
end


function create_kite_point_group(idx, points, wing, dynamics_type; simple=false)
    gamma = (-1 + 1/wing.n_groups + 2(idx-1)/wing.n_groups) * wing.gamma_tip
    le_pos = [wing.le_interp[i](gamma) for i in 1:3]
    chord = [wing.te_interp[i](gamma) for i in 1:3] .- le_pos
    y_airf = normalize([wing.le_interp[i](gamma-0.01) for i in 1:3] - le_pos)
    if simple
        top_point = wing.R_cad_body * (set.top_bridle_points[1] + wing.T_cad_body)
        gamma > 0 && (top_point .= top_point .* [1,-1,1])
        return KitePointGroup(idx, points, top_point, chord, y_airf, dynamics_type)
    else
        return KitePointGroup(idx, points, le_pos, chord, y_airf, dynamics_type)
    end
end

function create_kite_point(points, idx, set, wing; simple=false)
    points_per_group = length(set.bridle_fracs)
    kite_point_idx = count(p -> p.type == KITE, points) + 1
    group_idx = ceil(Int, kite_point_idx/points_per_group)
    gamma = (-1 + 1/wing.n_groups + 2(group_idx-1)/wing.n_groups) * wing.gamma_tip
    le_pos = [wing.le_interp[i](gamma) for i in 1:3]
    chord = [wing.te_interp[i](gamma) for i in 1:3] .- le_pos

    frac_idx = (kite_point_idx-1)%wing.n_groups+1
    if simple
        top_point = wing.R_cad_body * (set.top_bridle_points[1] + wing.T_cad_body)
        gamma > 0 && (top_point .= top_point .* [1,-1,1])
        pos = top_point .+ chord .* set.bridle_fracs[frac_idx]
    else
        pos = le_pos .+ chord .* set.bridle_fracs[frac_idx]
    end
    return Point(idx, pos, KITE)
end

function create_tether(idx, set, points, segments, tethers, attach_point, type, dynamics_type)
    l0 = set.l_tether / set.segments
    segment_idxs = Int16[]
    for i in 1:set.segments
        frac = i / set.segments
        pos = [(1-frac) * attach_point.pos_b[1], 
                (1-frac) * attach_point.pos_b[2],
                attach_point.pos_b[3] - i*l0]
        i_pnt = length(points) # last point idx
        i_seg = length(segments) # last segment idx
        if i == 1
            points = [points; Point(1+i_pnt, pos, dynamics_type)]
            segments = [segments; Segment(1+i_seg, (attach_point.idx, 1+i_pnt), type)]
        elseif i == set.segments
            points = [points; Point(1+i_pnt, pos, WINCH)]
            segments = [segments; Segment(1+i_seg, (i_pnt, 1+i_pnt), type)]
        else
            points = [points; Point(1+i_pnt, pos, dynamics_type)]
            segments = [segments; Segment(1+i_seg, (i_pnt, 1+i_pnt), type)]
        end
        push!(segment_idxs, 1+i_seg)
        i_pnt = length(points)
    end
    winch_point_idx = points[end].idx
    tethers = [tethers; Tether(idx, segment_idxs, winch_point_idx)]
    return points, segments, tethers, tethers[end].idx
end

function create_ram_point_system(set::Settings, wing::RamAirWing)
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
    
    bridle_top_left = [wing.R_cad_body * (set.top_bridle_points[i] + wing.T_cad_body) for i in eachindex(set.top_bridle_points)] # cad to kite frame
    bridle_top_right = [bridle_top_left[i] .* [1, -1, 1] for i in eachindex(set.top_bridle_points)]

    dynamics_type = set.quasi_static ? STATIC : DYNAMIC

    function create_bridle(bridle_top, gammas)
        i_pnt = length(points) # last point idx
        i_seg = length(segments) # last segment idx
        i_pul = length(pulleys) # last pulley idx

        i = 1
        for _ in gammas # 2 gammas per bridle system
            point_idxs = Int16[]
            for _ in set.bridle_fracs # 4 fracs
                points = [points; create_kite_point(points, i+i_pnt, set, wing)]
                i += 1
            end
            i_grp = 1 + length(groups)
            point_idxs = collect(points[end-3].idx:points[end].idx)
            groups = [groups; create_kite_point_group(i_grp, point_idxs, wing, DYNAMIC)]
        end

        points = [
            points
            Point(9+i_pnt, bridle_top[1], dynamics_type)
            Point(10+i_pnt, bridle_top[2], dynamics_type)
            Point(11+i_pnt, bridle_top[3], dynamics_type)
            Point(12+i_pnt, bridle_top[4], dynamics_type)

            Point(13+i_pnt, bridle_top[2] .+ [0, 0, -1], dynamics_type)

            Point(14+i_pnt, bridle_top[1] .+ [0, 0, -2], dynamics_type)
            Point(15+i_pnt, bridle_top[3] .+ [0, 0, -2], dynamics_type)

            Point(16+i_pnt, bridle_top[1] .+ [0, 0, -3], dynamics_type)
            Point(17+i_pnt, bridle_top[3] .+ [0, 0, -3], dynamics_type)
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
            Pulley(1+i_pul, (13+i_seg, 14+i_seg), dynamics_type)
            Pulley(2+i_pul, (16+i_seg, 17+i_seg), dynamics_type)
        ]
        push!(attach_points, points[end-1])
        push!(attach_points, points[end])
        return nothing
    end

    function create_tether(attach_point, type)
        l0 = set.l_tether / set.segments
        segment_idxs = Int16[]
        for i in 1:set.segments
            frac = i / set.segments
            pos = [(1-frac) * attach_point.pos_b[1], 
                    (1-frac) * attach_point.pos_b[2],
                    attach_point.pos_b[3] - i*l0]
            i_pnt = length(points) # last point idx
            i_seg = length(segments) # last segment idx
            if i == 1
                points = [points; Point(1+i_pnt, pos, dynamics_type)]
                segments = [segments; Segment(1+i_seg, (attach_point.idx, 1+i_pnt), type)]
            elseif i == set.segments
                points = [points; Point(1+i_pnt, pos, WINCH)]
                segments = [segments; Segment(1+i_seg, (i_pnt, 1+i_pnt), type)]
            else
                points = [points; Point(1+i_pnt, pos, dynamics_type)]
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

    winches = [winches; Winch(1, TorqueControlledMachine(set), [left_power_idx, right_power_idx])]
    winches = [winches; Winch(2, TorqueControlledMachine(set), [left_steering_idx])]
    winches = [winches; Winch(3, TorqueControlledMachine(set), [right_steering_idx])]

    return PointMassSystem(set.physical_model, points, groups, segments, pulleys, tethers, winches)
end

function create_simple_ram_point_system(set::Settings, wing::RamAirWing)
    (length(set.bridle_fracs) != 2) && throw(ArgumentError("Only 2 bridle fracs should be provided for the simple model."))

    points = Point[]
    groups = KitePointGroup[]
    segments = Segment[]
    pulleys = Pulley[]
    tethers = Tether[]
    winches = Winch[]

    dynamics_type = set.quasi_static ? STATIC : DYNAMIC

    [points = [points; create_kite_point(points, i, set, wing; simple=true)] for i in 1:4]

    groups = [
        groups
        create_kite_point_group(1, [1,2], wing, DYNAMIC; simple=true)
        create_kite_point_group(2, [3,4], wing, DYNAMIC; simple=true)
    ]

    points, segments, tethers, left_power_idx = create_tether(1, set, points, segments, tethers, points[1], POWER, dynamics_type)
    points, segments, tethers, right_power_idx = create_tether(2, set, points, segments, tethers, points[3], POWER, dynamics_type)
    points, segments, tethers, left_steering_idx = create_tether(3, set, points, segments, tethers, points[2], STEERING, dynamics_type)
    points, segments, tethers, right_steering_idx = create_tether(4, set, points, segments, tethers, points[4], STEERING, dynamics_type)

    winches = [winches; Winch(1, TorqueControlledMachine(set), [left_power_idx, right_power_idx])]
    winches = [winches; Winch(2, TorqueControlledMachine(set), [left_steering_idx])]
    winches = [winches; Winch(3, TorqueControlledMachine(set), [right_steering_idx])]

    return PointMassSystem(set.physical_model, points, groups, segments, pulleys, tethers, winches)
end


function init!(system::PointMassSystem, set::Settings, R_b_w)
    @unpack points, groups, segments, pulleys, tethers, winches = system

    for segment in segments
        (segment.type === BRIDLE) && (segment.diameter = 0.001set.bridle_tether_diameter)
        (segment.type === POWER) && (segment.diameter = 0.001set.power_tether_diameter)
        (segment.type === STEERING) && (segment.diameter = 0.001set.steering_tether_diameter)
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
        point.pos_w .= R_b_w * [point.pos_b[1], point.pos_b[2], point.pos_b[3] - min_z]
    end
    init_kite_pos = R_b_w * [0.0, 0.0, -min_z]
    return init_kite_pos
end

const MeasureFloat = Float32

@with_kw mutable struct Measurement
    set_values::MVector{3, MeasureFloat}    = [-50., -1., -1.]
    tether_length::MVector{3, MeasureFloat} = [51., 51., 51.]
    tether_vel::MVector{3, MeasureFloat}    = zeros(MeasureFloat, 3)
    tether_acc::MVector{3, MeasureFloat}    = zeros(MeasureFloat, 3)
    tether_force::MVector{3, MeasureFloat}  = [540., 3., 3.]
    "elevation and azimuth in spherical coordinate system with columns (left, right) and rows (elevation, azimuth)"
    sphere_pos::Matrix{MeasureFloat}            = deg2rad.([89.0 89.0; 1.0 -1.0])
    sphere_vel::Matrix{MeasureFloat}            = zeros(MeasureFloat, 2, 2)
    sphere_acc::Matrix{MeasureFloat}            = zeros(MeasureFloat, 2, 2)
    "positive azimuth wind direction in right-handed ENU frame relative to east / x-axis"
    wind_dir_gnd::MeasureFloat                  = zero(MeasureFloat)
end

function Base.getproperty(m::Measurement, val::Symbol)
    if val === :elevation
        sphere_pos = getfield(m, :sphere_pos)
        return 0.5(sphere_pos[1, 1] + sphere_pos[1, 2])
    elseif val === :azimuth
        sphere_pos = getfield(m, :sphere_pos)
        return 0.5(sphere_pos[2, 1] + sphere_pos[2, 2])
    else
        return getfield(m, val)
    end
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

