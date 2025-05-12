
function VortexStepMethod.RamAirWing(set::Settings; prn=true, kwargs...)
    obj_path = joinpath(dirname(get_data_path()), set.model)
    dat_path = joinpath(dirname(get_data_path()), set.foil_file)
    return RamAirWing(obj_path, dat_path; 
        mass=set.mass, crease_frac=set.crease_frac, n_groups=4, 
        align_to_principal=true, prn, kwargs...
    )
end

"""
    SegmentType `POWER` `STEERING` `BRIDLE`

Type of segment.

# Elements
- POWER: Belongs to a power line
- STEERING: Belongs to a steering line
- BRIDLE: Belongs to the bridle
"""
@enum SegmentType begin
    POWER
    STEERING
    BRIDLE
end

"""
    DynamicsType `DYNAMIC` `STATIC` `KITE` `WINCH`

Enumeration of the models that are attached to a point.

# Elements
- DYNAMIC: Belongs to a dynamic tether model
- STATIC: Belongs to a static tether model
- KITE: Rigid body
- WINCH: Winch
"""
@enum DynamicsType begin
    DYNAMIC
    STATIC
    KITE
    WINCH
end

"""
    mutable struct Point

A normal freely moving tether point.

$(TYPEDFIELDS)
"""
mutable struct Point
    idx::Int16
    pos_b::KVec3 # pos relative to kite COM in body frame
    pos_w::KVec3 # pos in world frame
    vel_w::KVec3 # vel in world frame
    type::DynamicsType
end
function Point(idx, pos_b, type, vel_w=zeros(KVec3))
    Point(idx, pos_b, copy(pos_b), vel_w, type)
end

"""
    struct KitePointGroup

Set of bridle lines that share the same twist angle and trailing edge angle.

$(TYPEDFIELDS)
"""
mutable struct KitePointGroup
    idx::Int16
    points::Vector{Int16}
    le_pos::KVec3 # point which the group rotates around under kite deformation
    chord::KVec3 # chord vector in body frame which the group rotates around under kite deformation
    y_airf::KVec3 # spanwise vector in local panel frame which the group rotates around under kite deformation
    type::DynamicsType
    moment_frac::SimFloat
    twist::SimFloat
    twist_vel::SimFloat
end
function KitePointGroup(idx, points, wing, gamma, type, moment_frac)
    le_pos = [wing.le_interp[i](gamma) for i in 1:3]
    chord = [wing.te_interp[i](gamma) for i in 1:3] .- le_pos
    y_airf = normalize([wing.le_interp[i](gamma-0.01) for i in 1:3] - le_pos)
    KitePointGroup(idx, points, le_pos, chord, y_airf, type, moment_frac, 0.0, 0.0)
end

"""
    mutable struct Segment

A segment from one point index to another point index.

$(TYPEDFIELDS)
"""
mutable struct Segment
    idx::Int16
    points::Tuple{Int16, Int16}
    type::SegmentType
    l0::SimFloat
    compression_frac::SimFloat
    diameter::SimFloat
end
function Segment(idx, points, type, l0=zero(SimFloat), compression_frac=0.1)
    Segment(idx, points, type, l0, compression_frac, zero(SimFloat))
end

"""
    mutable struct Pulley

A pulley described by two segments with the common point of the segments being the pulley.

$(TYPEDFIELDS)
"""
mutable struct Pulley
    idx::Int16
    segments::Tuple{Int16, Int16}
    type::DynamicsType
    sum_length::SimFloat
    length::SimFloat
    vel::SimFloat
    function Pulley(idx, segments, type)
        new(idx, segments, type, 0.0, 0.0, 0.0)
    end
end

"""
    struct Tether

A set of segments making a flexible tether. The winch point should only be part of one segment.

$(TYPEDFIELDS)
"""
struct Tether
    idx::Int16
    segments::Vector{Int16}
    winch_point::Int16
end

"""
    mutable struct Winch

A set of tethers or just one tether connected to a winch.

$(TYPEDFIELDS)
"""
mutable struct Winch
    idx::Int16
    model::AbstractWinchModel
    tethers::Vector{Int16}
    tether_length::SimFloat
    tether_vel::SimFloat
    function Winch(idx, model, tethers, tether_length, tether_vel=0.0)
        new(idx, model, tethers, tether_length, tether_vel)
    end
end

@with_kw struct Kite
    orient::KVec4 = zeros(KVec4)
    angular_vel::KVec3 = zeros(KVec3)
    pos::KVec3 = zeros(KVec3)
    vel::KVec3 = zeros(KVec3)
end

"""
    struct PointMassSystem

A discrete mass-spring-damper representation of a kite system, where point masses 
connected by elastic segments model the kite and tether dynamics:

- `points::Vector{Point}`: Point masses representing:
  - Kite attachment points 
  - Dynamic bridle/tether points
  - Fixed ground anchor points
- `groups::Vector{KitePointGroup}`: Collections of points that move together, 
    according to kite deformation (twist and trailing edge deflection)
- `segments::Vector{Segment}`: Spring-damper elements between points
- `pulleys::Vector{Pulley}`: Elements that redistribute line lengths
- `tethers::Vector{Tether}`: Groups of segments with a common unstretched length
- `winches::Vector{Winch}`: Ground-based winches that control the tether lengths

See also: [`Point`](@ref), [`Segment`](@ref), [`KitePointGroup`](@ref), [`Pulley`](@ref)
"""
struct PointMassSystem
    name::String
    points::Vector{Point}
    groups::Vector{KitePointGroup}
    segments::Vector{Segment}
    pulleys::Vector{Pulley}
    tethers::Vector{Tether}
    winches::Vector{Winch}
    kite::Kite
    function PointMassSystem(name, points, groups, segments, pulleys, tethers, winches, kite=Kite())
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
        return new(name, points, groups, segments, pulleys, tethers, winches, kite)
    end
end


function PointMassSystem(set::Settings, wing::RamAirWing)
    length(set.bridle_fracs) != 4 && throw(ArgumentError("4 bridle fracs should be provided for all models."))

    if set.physical_model == "ram"
        return create_ram_point_system(set, wing)
    elseif set.physical_model == "simple_ram"
        return create_simple_ram_point_system(set, wing)
    else
        throw(ArgumentError("Undefined physical model"))
    end
end

function calc_pos(wing::RamAirWing, gamma, frac)
    le_pos = [wing.le_interp[i](gamma) for i in 1:3]
    chord = [wing.te_interp[i](gamma) for i in 1:3] .- le_pos
    pos = le_pos .+ chord .* frac
    return pos
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
        i_grp = length(groups) # last group idx

        # ==================== CREATE DEFORMING KITE GROUPS ==================== #
        points = [
            points
            Point(1+i_pnt, calc_pos(wing, gammas[1], set.bridle_fracs[1]), KITE)
            Point(2+i_pnt, calc_pos(wing, gammas[1], set.bridle_fracs[3]), KITE)
            Point(3+i_pnt, calc_pos(wing, gammas[1], set.bridle_fracs[4]), KITE)

            Point(4+i_pnt, calc_pos(wing, gammas[2], set.bridle_fracs[1]), KITE)
            Point(5+i_pnt, calc_pos(wing, gammas[2], set.bridle_fracs[3]), KITE)
            Point(6+i_pnt, calc_pos(wing, gammas[2], set.bridle_fracs[4]), KITE)
        ]
        groups = [
            groups
            KitePointGroup(1+i_grp, [1+i_pnt, 2+i_pnt, 3+i_pnt], wing, gammas[1], DYNAMIC, set.bridle_fracs[2])
            KitePointGroup(2+i_grp, [4+i_pnt, 5+i_pnt, 6+i_pnt], wing, gammas[2], DYNAMIC, set.bridle_fracs[2])
        ]

        # ==================== CREATE PULLEY BRIDLE SYSTEM ==================== #
        points = [
            points
            Point(7+i_pnt, bridle_top[1], dynamics_type)
            Point(8+i_pnt, bridle_top[2], KITE)
            Point(9+i_pnt, bridle_top[3], dynamics_type)
            Point(10+i_pnt, bridle_top[4], dynamics_type)

            Point(11+i_pnt, bridle_top[2] .+ [0, 0, -1], dynamics_type)

            Point(12+i_pnt, bridle_top[1] .+ [0, 0, -2], dynamics_type)
            Point(13+i_pnt, bridle_top[3] .+ [0, 0, -2], dynamics_type)

            Point(14+i_pnt, bridle_top[1] .+ [0, 0, -4], dynamics_type)
            Point(15+i_pnt, bridle_top[3] .+ [0, 0, -4], dynamics_type)
        ]
        segments = [
            segments
            Segment(1+i_seg, (1+i_pnt, 7+i_pnt), BRIDLE)
            Segment(2+i_seg, (2+i_pnt, 9+i_pnt), BRIDLE)
            Segment(3+i_seg, (3+i_pnt, 10+i_pnt), BRIDLE)

            Segment(4+i_seg, (4+i_pnt, 7+i_pnt), BRIDLE)
            Segment(5+i_seg, (5+i_pnt, 9+i_pnt), BRIDLE)
            Segment(6+i_seg, (6+i_pnt, 10+i_pnt), BRIDLE)

            Segment(7+i_seg, (7+i_pnt, 12+i_pnt), BRIDLE, 2)
            Segment(8+i_seg, (8+i_pnt, 11+i_pnt), BRIDLE, 1)
            Segment(9+i_seg, (9+i_pnt, 13+i_pnt), BRIDLE, 2)
            Segment(10+i_seg, (10+i_pnt, 15+i_pnt), BRIDLE, 4)
            
            Segment(11+i_seg, (11+i_pnt, 12+i_pnt), BRIDLE, 1)
            Segment(12+i_seg, (11+i_pnt, 13+i_pnt), BRIDLE, 1)
            
            Segment(13+i_seg, (12+i_pnt, 14+i_pnt), BRIDLE, 2)
            Segment(14+i_seg, (13+i_pnt, 14+i_pnt), BRIDLE, 2)
            Segment(15+i_seg, (13+i_pnt, 15+i_pnt), BRIDLE, 2)
        ]
        pulleys = [
            pulleys
            Pulley(1+i_pul, (11+i_seg, 12+i_seg), dynamics_type)
            Pulley(2+i_pul, (14+i_seg, 15+i_seg), dynamics_type)
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
    create_bridle(bridle_top_left, gammas[[1,2]])
    create_bridle(bridle_top_right, gammas[[4,3]])

    left_power_idx = create_tether(attach_points[1], POWER)
    right_power_idx = create_tether(attach_points[3], POWER)
    left_steering_idx = create_tether(attach_points[2], STEERING)
    right_steering_idx = create_tether(attach_points[4], STEERING)

    winches = [winches; Winch(1, TorqueControlledMachine(set), [left_power_idx, right_power_idx], set.l_tether)]
    winches = [winches; Winch(2, TorqueControlledMachine(set), [left_steering_idx], set.l_tether)]
    winches = [winches; Winch(3, TorqueControlledMachine(set), [right_steering_idx], set.l_tether)]

    return PointMassSystem(set.physical_model, points, groups, segments, pulleys, tethers, winches)
end

function create_simple_ram_point_system(set::Settings, wing::RamAirWing)
    points = Point[]
    groups = KitePointGroup[]
    segments = Segment[]
    pulleys = Pulley[]
    tethers = Tether[]
    winches = Winch[]

    dynamics_type = set.quasi_static ? STATIC : DYNAMIC
    gammas = [-3/4, -1/4, 1/4, 3/4] * wing.gamma_tip
    
    bridle_top_left = [wing.R_cad_body * (set.top_bridle_points[i] + wing.T_cad_body) for i in eachindex(set.top_bridle_points)] # cad to kite frame
    bridle_top_right = [bridle_top_left[i] .* [1, -1, 1] for i in eachindex(set.top_bridle_points)]

    # ==================== CREATE DEFORMING KITE GROUPS ==================== #
    points = [
        points
        Point(1, calc_pos(wing, gammas[1], set.bridle_fracs[4]), KITE)
        Point(2, calc_pos(wing, gammas[2], set.bridle_fracs[4]), KITE)
        Point(3, calc_pos(wing, gammas[3], set.bridle_fracs[4]), KITE)
        Point(4, calc_pos(wing, gammas[4], set.bridle_fracs[4]), KITE)
    ]
    groups = [
        groups
        KitePointGroup(1, [1], wing, gammas[1], DYNAMIC, set.bridle_fracs[2])
        KitePointGroup(2, [2], wing, gammas[2], DYNAMIC, set.bridle_fracs[2])
        KitePointGroup(3, [3], wing, gammas[3], DYNAMIC, set.bridle_fracs[2])
        KitePointGroup(4, [4], wing, gammas[4], DYNAMIC, set.bridle_fracs[2])
    ]
    # ==================== CREATE PULLEY BRIDLE SYSTEM ==================== #
    points = [
        points
        Point(5, bridle_top_left[2], KITE)
        Point(6, bridle_top_left[4], dynamics_type)
        Point(7, bridle_top_right[2], KITE)
        Point(8, bridle_top_right[4], dynamics_type)
    ]

    segments = [
        segments
        Segment(1, (1, 6), BRIDLE)
        Segment(2, (2, 6), BRIDLE)
        Segment(3, (3, 8), BRIDLE)
        Segment(4, (4, 8), BRIDLE)
    ]

    points, segments, tethers, left_power_idx = create_tether(1, set, points, segments, tethers, points[5], POWER, dynamics_type)
    points, segments, tethers, right_power_idx = create_tether(2, set, points, segments, tethers, points[7], POWER, dynamics_type)
    points, segments, tethers, left_steering_idx = create_tether(3, set, points, segments, tethers, points[6], STEERING, dynamics_type)
    points, segments, tethers, right_steering_idx = create_tether(4, set, points, segments, tethers, points[8], STEERING, dynamics_type)

    winches = [winches; Winch(1, TorqueControlledMachine(set), [left_power_idx, right_power_idx], set.l_tether)]
    winches = [winches; Winch(2, TorqueControlledMachine(set), [left_steering_idx], set.l_tether)]
    winches = [winches; Winch(3, TorqueControlledMachine(set), [right_steering_idx], set.l_tether)]

    return PointMassSystem(set.physical_model, points, groups, segments, pulleys, tethers, winches)
end


function init!(system::PointMassSystem, set::Settings, R_b_w, Q_b_w)
    @unpack points, groups, segments, pulleys, tethers, winches, kite = system

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
        pulley.length = segment1.l0
        pulley.vel = 0.0
        @assert !(pulley.sum_length ≈ 0)
    end

    for winch in winches
        winch.tether_length = set.l_tether
        winch.tether_vel = 0.0
        @assert !(winch.tether_length ≈ 0)
    end

    first_moment_frac = groups[1].moment_frac
    for group in groups
        group.twist = 0.0
        group.twist_vel = 0.0
        @assert group.moment_frac ≈ first_moment_frac "All group.moment_frac must be the same."
    end

    min_z = Inf
    for point in points
        if point.pos_b[3] < min_z
            min_z = point.pos_b[3]
        end
    end
    for point in points
        point.pos_w .= R_b_w * [point.pos_b[1], point.pos_b[2], point.pos_b[3] - min_z]
        point.vel_w .= 0.0
    end
    kite.pos .= R_b_w * [0.0, 0.0, -min_z]
    kite.orient .= Q_b_w
    kite.vel .= 0.0
    kite.angular_vel .= 0.0
    return nothing
end

const MeasureFloat = Float32

@with_kw mutable struct Measurement
    set_values::MVector{3, MeasureFloat}    = [-50., -1., -1.]
    tether_length::MVector{3, MeasureFloat} = [51., 51., 51.]
    tether_vel::MVector{3, MeasureFloat}    = zeros(MeasureFloat, 3)
    tether_acc::MVector{3, MeasureFloat}    = zeros(MeasureFloat, 3)
    tether_force::MVector{3, MeasureFloat}  = [540., 3., 3.]
    "elevation and azimuth in spherical coordinate system with columns (left, right) and rows (elevation, azimuth)"
    sphere_pos::Matrix{MeasureFloat}            = deg2rad.([80.0 80.0; 1.0 -1.0])
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

function measure_to_q(measure::Measurement, R_cad_body=I(3))
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

