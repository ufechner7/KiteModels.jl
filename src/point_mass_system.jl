
# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

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
    DynamicsType `DYNAMIC` `QUASI_STATIC` `WING` `STATIC`

Enumeration of the models that are attached to a point.

# Elements
- DYNAMIC: Belongs to a dynamic tether model
- QUASI_STATIC: Belongs to a quasi static tether model
- WING: Connected to the rigid wing body
- STATIC: Does not change position
"""
@enum DynamicsType begin
    DYNAMIC
    QUASI_STATIC
    WING
    STATIC
end

"""
    mutable struct Point

A normal freely moving tether point.

$(TYPEDFIELDS)
"""
mutable struct Point
    idx::Int16
    wing_idx::Int16 # idx of wing used for initial orientation
    pos_b::KVec3 # pos relative to wing COM in body frame
    pos_w::KVec3 # pos in world frame
    vel_w::KVec3 # vel in world frame
    type::DynamicsType
end
function Point(idx, pos_b, type; vel_w=zeros(KVec3), wing_idx=1)
    Point(idx, wing_idx, pos_b, copy(pos_b), vel_w, type)
end

"""
    struct Group

Set of bridle lines that share the same twist angle and trailing edge angle.

$(TYPEDFIELDS)
"""
mutable struct Group
    idx::Int16
    point_idxs::Vector{Int16}
    le_pos::KVec3 # point which the group rotates around under wing deformation
    chord::KVec3 # chord vector in body frame which the group rotates around under wing deformation
    y_airf::KVec3 # spanwise vector in local panel frame which the group rotates around under wing deformation
    type::DynamicsType
    moment_frac::SimFloat
    twist::SimFloat
    twist_vel::SimFloat
end
function Group(idx, point_idxs, vsm_wing::RamAirWing, gamma, type, moment_frac)
    le_pos = [vsm_wing.le_interp[i](gamma) for i in 1:3]
    chord = [vsm_wing.te_interp[i](gamma) for i in 1:3] .- le_pos
    y_airf = normalize([vsm_wing.le_interp[i](gamma-0.01) for i in 1:3] - le_pos)
    Group(idx, point_idxs, le_pos, chord, y_airf, type, moment_frac, 0.0, 0.0)
end
function Group(idx, point_idxs, le_pos, chord, y_airf, type, moment_frac)
    Group(idx, point_idxs, le_pos, chord, y_airf, type, moment_frac, 0.0, 0.0)
end

"""
    mutable struct Segment

A segment from one point index to another point index.

$(TYPEDFIELDS)
"""
mutable struct Segment
    idx::Int16
    point_idxs::Tuple{Int16, Int16}
    type::SegmentType
    l0::SimFloat
    compression_frac::SimFloat
    diameter::SimFloat
end
function Segment(idx, point_idxs, type; l0=zero(SimFloat), compression_frac=0.1)
    Segment(idx, point_idxs, type, l0, compression_frac, zero(SimFloat))
end

"""
    mutable struct Pulley

A pulley described by two segments with the common point of the segments being the pulley.

$(TYPEDFIELDS)
"""
mutable struct Pulley
    idx::Int16
    segment_idxs::Tuple{Int16, Int16}
    type::DynamicsType
    sum_length::SimFloat
    length::SimFloat
    vel::SimFloat
    function Pulley(idx, segment_idxs, type)
        new(idx, segment_idxs, type, 0.0, 0.0, 0.0)
    end
end

"""
    struct Tether

A set of segments making a flexible tether. The winch point should only be part of one segment.

$(TYPEDFIELDS)
"""
struct Tether
    idx::Int16
    segment_idxs::Vector{Int16}
end

"""
    mutable struct Winch

A set of tethers or just one tether connected to a winch.

$(TYPEDFIELDS)
"""
mutable struct Winch
    idx::Int16
    model::AbstractWinchModel
    tether_idxs::Vector{Int16}
    tether_length::SimFloat
    tether_vel::SimFloat
    function Winch(idx, model, tether_idxs, tether_length, tether_vel=0.0)
        new(idx, model, tether_idxs, tether_length, tether_vel)
    end
end

@with_kw struct Wing
    idx::Int16
    group_idxs::Vector{Int16}
    orient::KVec4 = zeros(KVec4)
    angular_vel::KVec3 = zeros(KVec3)
    pos::KVec3 = zeros(KVec3)
    vel::KVec3 = zeros(KVec3)
end
function Wing(idx, group_idxs)
    Wing(; idx, group_idxs)
end

"""
    struct SystemStructure

A discrete mass-spring-damper representation of a kite system, where point masses 
connected by elastic segments model the kite and tether dynamics:

- `points::Vector{Point}`: Point masses representing:
  - Wing attachment points 
  - Dynamic bridle/tether points
  - Fixed ground anchor points
- `groups::Vector{Group}`: Collections of points that move together, 
    according to wing deformation (twist and trailing edge deflection)
- `segments::Vector{Segment}`: Spring-damper elements between points
- `pulleys::Vector{Pulley}`: Elements that redistribute line lengths
- `tethers::Vector{Tether}`: Groups of segments with a common unstretched length
- `winches::Vector{Winch}`: Ground-based winches that control the tether lengths

See also: [`Point`](@ref), [`Segment`](@ref), [`Group`](@ref), [`Pulley`](@ref)
"""
struct SystemStructure
    name::String
    points::Vector{Point}
    groups::Vector{Group}
    segments::Vector{Segment}
    pulleys::Vector{Pulley}
    tethers::Vector{Tether}
    winches::Vector{Winch}
    wings::Vector{Wing}
    y::Array{Float64, 2}
    x::Array{Float64, 2}
    jac::Array{Float64, 3}
    function SystemStructure(name; 
            points=Point[], 
            groups=Group[], 
            segments=Segment[], 
            pulleys=Pulley[], 
            tethers=Tether[], 
            winches=Winch[], 
            wings=Wing[]
        )
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
        if length(wings) > 0
            ny = 3+length(wings[1].group_idxs)+3
            nx = 3+3+length(wings[1].group_idxs)
        else
            ny = 0
            nx = 0
        end
        y = zeros(length(wings), ny)
        x = zeros(length(wings), nx)
        jac = zeros(length(wings), nx, ny)
        return new(name, points, groups, segments, pulleys, tethers, winches, wings, y, x, jac)
    end
end


function SystemStructure(set::Settings, wing::RamAirWing)
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

function create_tether(tether_idx, set, points, segments, tethers, attach_point, type, dynamics_type, z=[0,0,1])
    winch_pos = find_axis_point(attach_point.pos_b, set.l_tether, z)
    dir = winch_pos - attach_point.pos_b
    segment_idxs = Int16[]
    for i in 1:set.segments
        frac = i / set.segments
        pos = attach_point.pos_b + frac * dir
        point_idx = length(points)+1 # last point idx
        segment_idx = length(segments)+1 # last segment idx
        if i == 1
            points = [points; Point(point_idx, pos, dynamics_type)]
            segments = [segments; Segment(segment_idx, (attach_point.idx, point_idx), type)]
        elseif i == set.segments
            points = [points; Point(point_idx, pos, STATIC)]
            segments = [segments; Segment(segment_idx, (point_idx-1, point_idx), type)]
        else
            points = [points; Point(point_idx, pos, dynamics_type)]
            segments = [segments; Segment(segment_idx, (point_idx-1, point_idx), type)]
        end
        push!(segment_idxs, segment_idx)
    end
    tethers = [tethers; Tether(tether_idx, segment_idxs)]
    return points, segments, tethers, tethers[end].idx
end

function cad_to_body_frame(wing::RamAirWing, pos)
    return wing.R_cad_body * (pos + wing.T_cad_body)
end

# Find the point on the z-axis with distance l from P in the negative direction
# TODO: rename P to pos
function find_axis_point(P, l, v=[0,0,1])
    # Compute dot product v · P
    v ⋅ P = v[1] * P[1] + v[2] * P[2] + v[3] * P[3]
    # Compute discriminant
    D = (v ⋅ P)^2 - norm(v)^2 * (norm(P)^2 - l^2)
    D < 0 && error("No real solution: l is too small or parameters invalid")
    # Solve quadratic for t, choose solution for negative direction
    t = (v ⋅ P - √D) / norm(v)^2
    # Compute point Q = t * v
    return [t * v[1], t * v[2], t * v[3]]
end

function create_ram_point_system(set::Settings, vsm_wing::RamAirWing)
    points = Point[]
    groups = Group[]
    segments = Segment[]
    pulleys = Pulley[]
    tethers = Tether[]
    winches = Winch[]
    wings = Wing[]

    attach_points = Point[]
    
    bridle_top_left = [cad_to_body_frame(vsm_wing, set.top_bridle_points[i]) for i in eachindex(set.top_bridle_points)]
    bridle_top_right = [bridle_top_left[i] .* [1, -1, 1] for i in eachindex(set.top_bridle_points)]

    dynamics_type = set.quasi_static ? QUASI_STATIC : DYNAMIC
    z = vsm_wing.R_cad_body[:,3]

    function create_bridle(bridle_top, gammas)
        i_pnt = length(points) # last point idx
        i_seg = length(segments) # last segment idx
        i_pul = length(pulleys) # last pulley idx
        i_grp = length(groups) # last group idx

        # ==================== CREATE DEFORMING WING GROUPS ==================== #
        points = [
            points
            Point(1+i_pnt, calc_pos(vsm_wing, gammas[1], set.bridle_fracs[1]), WING)
            Point(2+i_pnt, calc_pos(vsm_wing, gammas[1], set.bridle_fracs[3]), WING)
            Point(3+i_pnt, calc_pos(vsm_wing, gammas[1], set.bridle_fracs[4]), WING)

            Point(4+i_pnt, calc_pos(vsm_wing, gammas[2], set.bridle_fracs[1]), WING)
            Point(5+i_pnt, calc_pos(vsm_wing, gammas[2], set.bridle_fracs[3]), WING)
            Point(6+i_pnt, calc_pos(vsm_wing, gammas[2], set.bridle_fracs[4]), WING)
        ]
        groups = [
            groups
            Group(1+i_grp, [1+i_pnt, 2+i_pnt, 3+i_pnt], vsm_wing, gammas[1], DYNAMIC, set.bridle_fracs[2])
            Group(2+i_grp, [4+i_pnt, 5+i_pnt, 6+i_pnt], vsm_wing, gammas[2], DYNAMIC, set.bridle_fracs[2])
        ]

        # ==================== CREATE PULLEY BRIDLE SYSTEM ==================== #
        # TODO: add initial rotation around y-axis
        points = [
            points
            Point(7+i_pnt, bridle_top[1], dynamics_type)
            Point(8+i_pnt, bridle_top[2], WING)
            Point(9+i_pnt, bridle_top[3], dynamics_type)
            Point(10+i_pnt, bridle_top[4], dynamics_type)

            Point(11+i_pnt, bridle_top[2] .+ -1z, dynamics_type)

            Point(12+i_pnt, bridle_top[1] .+ -2z, dynamics_type)
            Point(13+i_pnt, bridle_top[3] .+ -2z, dynamics_type)

            Point(14+i_pnt, bridle_top[1] .+ -4z, dynamics_type)
            Point(15+i_pnt, bridle_top[3] .+ -4z, dynamics_type)
        ]
        segments = [
            segments
            Segment(1+i_seg, (1+i_pnt, 7+i_pnt), BRIDLE)
            Segment(2+i_seg, (2+i_pnt, 9+i_pnt), BRIDLE)
            Segment(3+i_seg, (3+i_pnt, 10+i_pnt), BRIDLE)

            Segment(4+i_seg, (4+i_pnt, 7+i_pnt), BRIDLE)
            Segment(5+i_seg, (5+i_pnt, 9+i_pnt), BRIDLE)
            Segment(6+i_seg, (6+i_pnt, 10+i_pnt), BRIDLE)

            Segment(7+i_seg, (7+i_pnt, 12+i_pnt), BRIDLE; l0=2)
            Segment(8+i_seg, (8+i_pnt, 11+i_pnt), BRIDLE; l0=1)
            Segment(9+i_seg, (9+i_pnt, 13+i_pnt), BRIDLE; l0=2)
            Segment(10+i_seg, (10+i_pnt, 15+i_pnt), BRIDLE; l0=4)
            
            Segment(11+i_seg, (11+i_pnt, 12+i_pnt), BRIDLE; l0=1)
            Segment(12+i_seg, (11+i_pnt, 13+i_pnt), BRIDLE; l0=1)
            
            Segment(13+i_seg, (12+i_pnt, 14+i_pnt), BRIDLE; l0=2)
            Segment(14+i_seg, (13+i_pnt, 14+i_pnt), BRIDLE; l0=2)
            Segment(15+i_seg, (13+i_pnt, 15+i_pnt), BRIDLE; l0=2)
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

    gammas = [-3/4, -1/4, 1/4, 3/4] * vsm_wing.gamma_tip
    create_bridle(bridle_top_left, gammas[[1,2]])
    create_bridle(bridle_top_right, gammas[[3,4]])

    points, segments, tethers, left_power_idx = create_tether(1, set, points, segments, tethers, attach_points[1], POWER, dynamics_type, z)
    points, segments, tethers, right_power_idx = create_tether(2, set, points, segments, tethers, attach_points[3], POWER, dynamics_type, z)
    points, segments, tethers, left_steering_idx = create_tether(3, set, points, segments, tethers, attach_points[2], STEERING, dynamics_type, z)
    points, segments, tethers, right_steering_idx = create_tether(4, set, points, segments, tethers, attach_points[4], STEERING, dynamics_type, z)

    winches = [winches; Winch(1, TorqueControlledMachine(set), [left_power_idx, right_power_idx], set.l_tether)]
    winches = [winches; Winch(2, TorqueControlledMachine(set), [left_steering_idx], set.l_tether)]
    winches = [winches; Winch(3, TorqueControlledMachine(set), [right_steering_idx], set.l_tether)]

    wings = [Wing(1, [1,2,3,4])]
    
    return SystemStructure(set.physical_model; points, groups, segments, pulleys, tethers, winches, wings)
end

function create_simple_ram_point_system(set::Settings, wing::RamAirWing)
    points = Point[]
    groups = Group[]
    segments = Segment[]
    pulleys = Pulley[]
    tethers = Tether[]
    winches = Winch[]

    dynamics_type = set.quasi_static ? QUASI_STATIC : DYNAMIC
    gammas = [-3/4, -1/4, 1/4, 3/4] * wing.gamma_tip
    
    bridle_top_left = [wing.R_cad_body * (set.top_bridle_points[i] + wing.T_cad_body) for i in eachindex(set.top_bridle_points)] # cad to kite frame
    bridle_top_right = [bridle_top_left[i] .* [1, -1, 1] for i in eachindex(set.top_bridle_points)]

    # ==================== CREATE DEFORMING WING GROUPS ==================== #
    points = [
        points
        Point(1, calc_pos(wing, gammas[1], set.bridle_fracs[4]), WING)
        Point(2, calc_pos(wing, gammas[2], set.bridle_fracs[4]), WING)
        Point(3, calc_pos(wing, gammas[3], set.bridle_fracs[4]), WING)
        Point(4, calc_pos(wing, gammas[4], set.bridle_fracs[4]), WING)
    ]
    groups = [
        groups
        Group(1, [1], wing, gammas[1], DYNAMIC, set.bridle_fracs[2])
        Group(2, [2], wing, gammas[2], DYNAMIC, set.bridle_fracs[2])
        Group(3, [3], wing, gammas[3], DYNAMIC, set.bridle_fracs[2])
        Group(4, [4], wing, gammas[4], DYNAMIC, set.bridle_fracs[2])
    ]
    # ==================== CREATE PULLEY BRIDLE SYSTEM ==================== #
    points = [
        points
        Point(5, bridle_top_left[2], WING)
        Point(6, bridle_top_left[4], dynamics_type)
        Point(7, bridle_top_right[2], WING)
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

    return SystemStructure(set.physical_model, points, groups, segments, pulleys, tethers, winches)
end


function init!(system::SystemStructure, set::Settings, R_b_w, Q_b_w)
    @unpack points, groups, segments, pulleys, tethers, winches, wings = system

    for segment in segments
        (segment.type === BRIDLE) && (segment.diameter = 0.001set.bridle_tether_diameter)
        (segment.type === POWER) && (segment.diameter = 0.001set.power_tether_diameter)
        (segment.type === STEERING) && (segment.diameter = 0.001set.steering_tether_diameter)
        (segment.l0 ≈ 0) && (segment.l0 = norm(points[segment.point_idxs[1]].pos_b - points[segment.point_idxs[2]].pos_b) * 0.9999)
        @assert (0 < segment.diameter < 1)
        @assert (segment.l0 > 0)
    end

    for pulley in pulleys
        segment1, segment2 = segments[pulley.segment_idxs[1]], segments[pulley.segment_idxs[2]]
        pulley.sum_length = segment1.l0 + segment2.l0
        pulley.length = segment1.l0
        pulley.vel = 0.0
        @assert !(pulley.sum_length ≈ 0)
    end

    for winch in winches
        winch.tether_length ≈ 0.0 && (winch.tether_length = set.l_tether)
        winch.tether_vel = 0.0
        @assert !(winch.tether_length ≈ 0)
    end

    (length(groups) > 0) && (first_moment_frac = groups[1].moment_frac)
    for group in groups
        group.twist = 0.0
        group.twist_vel = 0.0
        @assert group.moment_frac ≈ first_moment_frac "All group.moment_frac must be the same."
    end

    min_point = fill(Inf, 3)
    for point in points
        if point.pos_b[3] < min_point[3]
            min_point .= point.pos_b
        end
    end
    for point in points
        if point.wing_idx == 0
            R = I(3)
        else
            R = R_b_w[point.wing_idx, :, :]
        end
        point.pos_w .= R * (point.pos_b .- min_point)
        point.vel_w .= 0.0
    end
    for wing in wings
        wing.pos .= R_b_w[wing.idx, :, :] * -min_point
        wing.orient .= Q_b_w[wing.idx, :]
        wing.vel .= 0.0
        wing.angular_vel .= 0.0
    end
    return nothing
end
