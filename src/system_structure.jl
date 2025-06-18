
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
struct Point
    idx::Int16
    transform_idx::Int16 # idx of wing used for initial orientation
    wing_idx::Int16
    pos_cad::KVec3
    pos_b::KVec3 # pos relative to wing COM in body frame
    pos_w::KVec3 # pos in world frame
    vel_w::KVec3 # vel in world frame
    type::DynamicsType
end

"""
    Point(idx, pos_cad, type; wing_idx=1, vel_w=zeros(KVec3), transform_idx=1)

Constructs a Point object.

# Arguments
- `idx::Int16`: Unique identifier for the point.
- `pos_cad::KVec3`: Position of the point in the CAD frame.
- `type::DynamicsType`: Dynamics type of the point (STATIC, DYNAMIC, etc.).

# Keyword Arguments
- `wing_idx::Int16=1`: Index of the wing this point is attached to.
- `vel_w::KVec3=zeros(KVec3)`: Initial velocity of the point in world frame.
- `transform_idx::Int16=1`: Index of the transform used for initial positioning.

# Returns
- `Point`: A new Point object.

# Example
To create a Point:
```julia
    point = Point(1, [1.0, 2.0, 3.0], DYNAMIC; wing_idx=1)
```
"""
function Point(idx, pos_cad, type; wing_idx=1, vel_w=zeros(KVec3), transform_idx=1)
    Point(idx, transform_idx, wing_idx, pos_cad, zeros(KVec3), zeros(KVec3), vel_w, type)
end

"""
    struct Group

Set of bridle lines that share the same twist angle and trailing edge angle.

$(TYPEDFIELDS)
"""
mutable struct Group
    const idx::Int16
    const point_idxs::Vector{Int16}
    const le_pos::KVec3 # point which the group rotates around under wing deformation
    const chord::KVec3 # chord vector in body frame which the group rotates around under wing deformation
    const y_airf::KVec3 # spanwise vector in local panel frame which the group rotates around under wing deformation
    const type::DynamicsType
    moment_frac::SimFloat
    twist::SimFloat
    twist_vel::SimFloat
end

"""
    Group(idx, point_idxs, vsm_wing::RamAirWing, gamma, type, moment_frac)

Constructs a Group object representing a collection of points on a kite body that share 
a common twist deformation.

A Group models the local deformation of a kite wing section through twist dynamics. 
All points within a group undergo the same twist rotation about the chord vector, 
which is determined by a torque balance between aerodynamic moments and tether forces.

The governing equation is:
```julia
  twist_α = (group_aero_moment + group_tether_moment) / inertia
```

# Arguments
- `idx::Int16`: Unique identifier for the group
- `point_idxs::Vector{Int16}`: Indices of points that move together with this group's twist
- `vsm_wing::RamAirWing`: Wing geometry object used to extract local chord and spanwise vectors
- `gamma`: Spanwise parameter (typically -1 to 1) defining the group's location along the wing
- `type::DynamicsType`: Dynamics type (DYNAMIC for time-varying twist, QUASI_STATIC for equilibrium)
- `moment_frac::SimFloat`: Fraction of total wing moment applied to this group (typically distributed equally)

# Returns
- `Group`: A new Group object with twist dynamics capability

# Example
Create a group at mid-span with quarter of the wing moment:
```julia
  group = Group(1, [1, 2, 3], vsm_wing, 0.0, DYNAMIC, 0.25)
```
"""
function Group(idx, point_idxs, vsm_wing::RamAirWing, gamma, type, moment_frac)
    le_pos = [vsm_wing.le_interp[i](gamma) for i in 1:3]
    chord = [vsm_wing.te_interp[i](gamma) for i in 1:3] .- le_pos
    y_airf = normalize([vsm_wing.le_interp[i](gamma-0.01) for i in 1:3] - le_pos)
    Group(idx, point_idxs, le_pos, chord, y_airf, type, moment_frac, 0.0, 0.0)
end

"""
    Group(idx, point_idxs, le_pos, chord, y_airf, type, moment_frac)

Constructs a Group object representing a collection of points on a kite body that share 
a common twist deformation.

A Group models the local deformation of a kite wing section through twist dynamics. 
All points within a group undergo the same twist rotation about the chord vector, 
which is determined by a torque balance between aerodynamic moments and tether forces.

The governing equation is:
```julia
  twist_α = (group_aero_moment + group_tether_moment) / inertia
```

# Arguments
- `idx::Int16`: Unique identifier for the group
- `point_idxs::Vector{Int16}`: Indices of points that twist together with this group
- `le_pos::KVec3`: Leading edge position serving as the rotation center in body frame
- `chord::KVec3`: Chord vector defining the twist axis direction in body frame
- `y_airf::KVec3`: Spanwise vector in local airfoil frame for coordinate system definition
- `type::DynamicsType`: DYNAMIC for time-varying twist, QUASI_STATIC for equilibrium twist
- `moment_frac::SimFloat`: Fraction of total wing moment applied to this group

# Returns
- `Group`: A new Group object for modeling local wing deformation

# Example
Create a group with explicit geometry for a rectangular wing section:
```julia
  le_pos = [1.0, 2.0, 3.0]      # Leading edge position
  chord  = [1.0, 0.0, 0.0]       # Chord vector pointing downstream  
  y_airf = [0.0, 1.0, 0.0]      # Spanwise direction
  group = Group(1, [1, 2, 3], le_pos, chord, y_airf, DYNAMIC, 0.25)
```
"""
function Group(idx, point_idxs, le_pos, chord, y_airf, type, moment_frac)
    Group(idx, point_idxs, le_pos, chord, y_airf, type, moment_frac, 0.0, 0.0)
end

"""
    mutable struct Segment

A segment from one point index to another point index.

$(TYPEDFIELDS)
"""
mutable struct Segment
    const idx::Int16
    const point_idxs::Tuple{Int16, Int16}
    const type::SegmentType
    l0::SimFloat
    compression_frac::SimFloat
    diameter::SimFloat
end

"""
    Segment(idx, point_idxs, type; l0=zero(SimFloat), compression_frac=0.1)

Constructs a Segment object. If l0 is not provided, it is calculated as the distance between the two segment points.

# Arguments
- `idx::Int16`: Unique identifier for the segment.
- `point_idxs::Tuple{Int16, Int16}`: Tuple containing the indices of the two points connected by this segment.
- `type::SegmentType`: Type of the segment (POWER, STEERING, BRIDLE).

# Keyword Arguments
- `l0::SimFloat=zero(SimFloat)`: Unstretched length of the segment.
- `compression_frac::SimFloat=0.1`: Compression fraction of the segment.

# Returns
- `Segment`: A new Segment object.

# Example
To create a Segment:
```julia
    segment = Segment(1, (1, 2), BRIDLE; l0=10.0)
```
"""
function Segment(idx, point_idxs, type; l0=zero(SimFloat), compression_frac=0.1)
    Segment(idx, point_idxs, type, l0, compression_frac, zero(SimFloat))
end

"""
    mutable struct Pulley

A pulley described by two segments with the common point of the segments being the pulley.

$(TYPEDFIELDS)
"""
mutable struct Pulley
    const idx::Int16
    const segment_idxs::Tuple{Int16, Int16}
    const type::DynamicsType
    sum_length::SimFloat
    length::SimFloat
    vel::SimFloat
end
"""
    Pulley(idx, segment_idxs, type)

Constructs a Pulley object.

# Arguments
- `idx::Int16`: Unique identifier for the pulley.
- `segment_idxs::Tuple{Int16, Int16}`: Tuple containing the indices of the two segments connected by this pulley.
- `type::DynamicsType`: Dynamics type of the pulley.

# Returns
- `Pulley`: A new Pulley object.

# Example
To create a Pulley:
```julia
    pulley = Pulley(1, (1, 2), DYNAMIC)
```
"""
function Pulley(idx, segment_idxs, type)
    return Pulley(idx, segment_idxs, type, 0.0, 0.0, 0.0)
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
    const idx::Int16
    const model::AbstractWinchModel
    const tether_idxs::Vector{Int16}
    tether_length::SimFloat
    tether_vel::SimFloat
end

"""
    Winch(idx, model, tether_idxs, tether_length; tether_vel=0.0)

Constructs a Winch object.

# Arguments
- `idx::Int16`: Unique identifier for the winch.
- `model::AbstractWinchModel`: The winch model.
- `tether_idxs::Vector{Int16}`: Vector containing the indices of the tethers connected to this winch.
- `tether_length::SimFloat`: Initial tether length.

# Keyword Arguments
- `tether_vel::SimFloat=0.0`: Initial tether velocity.

# Returns
- `Winch`: A new Winch object.

# Example
To create a Winch:
```julia
    winch = Winch(1, TorqueControlledMachine(set), [1, 2], 100.0)
```
"""
function Winch(idx, model, tether_idxs, tether_length; tether_vel=0.0)
    return Winch(idx, model, tether_idxs, tether_length, tether_vel)
end

"""
    struct Wing

A rigid wing body that can have multiple groups of points attached to it.

# Fields
- `idx::Int16`: Unique identifier for the wing
- `group_idxs::Vector{Int16}`: Indices of groups attached to this wing
- `transform_idx::Int16`: Transform used for initial positioning and orientation
- `R_b_c::Matrix{SimFloat}`: Rotation matrix from body frame to CAD frame
- `angular_vel::KVec3`: Angular velocity of the wing in world frame
- `pos_w::KVec3`: Position of wing center of mass in world frame
- `pos_cad::KVec3`: Position of wing center of mass in CAD frame
- `vel_w::KVec3`: Velocity of wing center of mass in world frame

The wing provides a rigid body reference frame for attached points and groups.
Points with `type == WING` move rigidly with the wing body according to the
wing's orientation matrix `R_b_c` and position `pos_w`.

# Extended help
The wing's orientation can be accessed as a quaternion through the `orient` property:
```julia
wing = Wing(1, [1,2], I(3), zeros(3))
quat = wing.orient  # Returns quaternion representation of R_b_c
```
"""
struct Wing
    idx::Int16
    group_idxs::Vector{Int16}
    transform_idx::Int16
    R_b_c::Matrix{SimFloat}
    angular_vel::KVec3
    pos_w::KVec3
    pos_cad::KVec3
    vel_w::KVec3
end
function Base.getproperty(wing::Wing, sym::Symbol)
    if sym == :orient
        return rotation_matrix_to_quaternion(wing.R_b_c)
    else
        return getfield(wing, sym)
    end
end

"""
    Wing(idx, group_idxs, R_b_c, pos_cad; transform_idx=1, angular_vel=zeros(KVec3), 
         pos_w=zeros(KVec3), vel_w=zeros(KVec3))

Constructs a Wing object representing a rigid body that serves as a reference frame for attached points and groups.

A Wing provides a rigid body coordinate system for kite components. Points with `type == WING` move rigidly 
with the wing body according to the wing's orientation matrix and position. Groups attached to the wing 
undergo local deformation (twist) relative to the rigid wing body frame.

# Arguments
- `idx::Int16`: Unique identifier for the wing
- `group_idxs::Vector{Int16}`: Indices of groups attached to this wing that can deform relative to the body
- `R_b_c::Matrix{SimFloat}`: Rotation matrix from body frame to CAD frame (3×3 orthogonal matrix)
- `pos_cad::KVec3`: Position of wing center of mass in CAD frame

# Keyword Arguments
- `transform_idx::Int16=1`: Transform used for initial positioning and orientation
- `angular_vel::KVec3=zeros(KVec3)`: Initial angular velocity of the wing in world frame
- `pos_w::KVec3=zeros(KVec3)`: Initial position of wing center of mass in world frame
- `vel_w::KVec3=zeros(KVec3)`: Initial velocity of wing center of mass in world frame

# Special Properties
The wing orientation can be accessed as a quaternion:
```julia
  quat = wing.orient  # Returns quaternion representation of R_b_c
```

# Returns
- `Wing`: A new Wing object providing a rigid body reference frame

# Example
Create a wing with identity orientation and two attached groups:
```julia
  R_b_c = I(3) # identity matrix
  pos_cad = [0.0, 0.0, 0.0]
  wing = Wing(1, [1, 2], R_b_c, pos_cad)
```
"""
function Wing(idx, group_idxs, R_b_c, pos_cad; transform_idx=1, angular_vel=zeros(KVec3), 
        pos_w=zeros(KVec3), vel_w=zeros(KVec3))
    return Wing(idx, group_idxs, transform_idx, R_b_c, angular_vel, pos_w, pos_cad, vel_w)
end

"""
    mutable struct Transform

Describes the spatial transformation (position and orientation) of system components
relative to a base reference point.

$(TYPEDFIELDS)
"""
mutable struct Transform
    const idx::Int16
    const wing_idx::Union{Int16, Nothing}
    const rot_point_idx::Union{Int16, Nothing}
    const base_point_idx::Union{Int16, Nothing}
    const base_transform_idx::Union{Int16, Nothing}
    elevation::SimFloat # The elevation of the rotating point or kite as seen from the base point
    azimuth::SimFloat # The azimuth of the rotating point or kite as seen from the base point
    heading::SimFloat
    base_pos::Union{KVec3, Nothing}
end

"""
    Transform(idx, elevation, azimuth, heading; 
        base_point_idx=nothing, base_pos=nothing, base_transform_idx=nothing, 
        wing_idx=nothing, rot_point_idx=nothing)

Constructs a Transform object that describes an orientation transformation from a base reference to a rotating object.

A Transform defines the spatial relationship between a reference (base) and a target object (wing or point) 
using spherical coordinates. The transformation is defined by elevation, azimuth, and heading angles that position 
and orient the target relative to the base.

# Base Reference
The base of the transformation can be defined in two ways:
- **Fixed position**: Using `base_pos` (fixed position vector) and `base_point_idx` (reference point index)
- **Chained transform**: Using `base_transform_idx` (index of another transform whose position becomes the base)

# Target Object
The target of the transformation can be either:
- A wing (specified by `wing_idx`): The entire wing will be positioned and oriented
- A point (specified by `rot_point_idx`): A single point will be positioned

# Arguments
- `idx::Int16`: Unique identifier for the transform
- `elevation::SimFloat`: Elevation angle of the target as seen from base (radians)
- `azimuth::SimFloat`: Azimuth angle of the target as seen from base (radians)
- `heading::SimFloat`: Rotation angle around the base-to-target vector (radians)

# Keyword Arguments
- **Base Reference (mutually exclusive)**
- `base_pos::AbstractVector=nothing`: Fixed position offset for the base point in world coordinates
- `base_point_idx::Int16=nothing`: Index of the reference point (required when using `base_pos`)
- `base_transform_idx::Int=nothing`: Index of another transform to use as base position

- **Target Object (mutually exclusive)**
- `wing_idx::Union{Int16, Nothing}=nothing`: Index of wing to be rotated to (elevation, azimuth)
- `rot_point_idx::Union{Int16, Nothing}=nothing`: Index of point to be rotated to (elevation, azimuth)

# Returns
- `Transform`: A new Transform object defining the spatial relationship

# Examples
Position a wing at 45° elevation and 30° azimuth from a fixed position:
  transform = Transform(1, deg2rad(45), deg2rad(30), 0.0; 
                       base_pos=[0,0,0], base_point_idx=1, wing_idx=1)

Position a point relative to another transform:
  transform = Transform(2, deg2rad(60), 0.0, 0.0; 
                       base_transform_idx=1, rot_point_idx=5)

Create a chained wing positioning (wing relative to another wing's position):
  transform = Transform(3, deg2rad(30), deg2rad(45), deg2rad(10); 
                       base_transform_idx=1, wing_idx=2)
"""
function Transform(idx, elevation, azimuth, heading;
        base_point_idx=nothing, base_pos=nothing, base_transform_idx=nothing,
        wing_idx=nothing, rot_point_idx=nothing)
    (isnothing(wing_idx) == isnothing(rot_point_idx)) && error("Either provide a wing_idx or a rot_point_idx, not both or none.")
    (isnothing(base_pos) == isnothing(base_transform_idx)) && error("Either provide the base_pos or the base_transform_idx, not both or none.")
    (isnothing(base_pos) !== isnothing(base_point_idx)) && error("When providing a base_pos, also provide a base_point_idx.")
    Transform(idx, wing_idx, rot_point_idx, base_point_idx, base_transform_idx, elevation, azimuth, heading, base_pos)
end
function Transform(idx, set, base_point_idx; kwargs...)
    Transform(idx, set.elevations[idx], set.azimuths[idx], set.headings[idx], base_point_idx; kwargs...)
end

function get_rot_pos(transform::Transform, wings, points)
    if !isnothing(transform.wing_idx)
        return wings[transform.wing_idx].pos_w
    elseif !isnothing(transform.rot_point_idx)
        return points[transform.rot_point_idx].pos_w
    end
end

function get_base_pos(transform::Transform, wings, points)
    curr_base_pos = points[transform.base_point_idx].pos_cad
    if !isnothing(transform.base_pos)
        return transform.base_pos, curr_base_pos
    elseif !isnothing(transform.base_transform_idx)
        base_transform = transforms[transform.base_transform_idx]
        return get_rot_pos(base_transform, wings, points), curr_base_pos
    end
end

"""
    struct SystemStructure

A discrete mass-spring-damper representation of a kite sys_struct, where point masses 
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
    set::Settings
    points::Vector{Point}
    groups::Vector{Group}
    segments::Vector{Segment}
    pulleys::Vector{Pulley}
    tethers::Vector{Tether}
    winches::Vector{Winch}
    wings::Vector{Wing}
    transforms::Vector{Transform}
    y::Array{Float64, 2}
    x::Array{Float64, 2}
    jac::Array{Float64, 3}
end

"""
    SystemStructure(name, set; points=Point[], groups=Group[], segments=Segment[], 
                   pulleys=Pulley[], tethers=Tether[], winches=Winch[], 
                   wings=Wing[], transforms=Transform[])

Constructs a SystemStructure object representing a complete kite system using a discrete mass-spring-damper model.

A SystemStructure defines the physical topology and dynamics of an airborne wind energy system, where:
- Point masses represent wing attachment points, bridle connections, and tether nodes
- Spring-damper segments connect points to model tether/bridle elasticity
- Groups enable local wing deformation through twist dynamics
- Pulleys redistribute line tensions and lengths
- Winches provide tether actuation and control

# Physical Model
The system uses a lumped-parameter approach where:
- **Points**: Discrete masses with 3DOF translational dynamics
- **Segments**: Elastic connections with spring-damper behavior
- **Groups**: Collections of points sharing twist deformation
- **Wings**: Rigid bodies providing reference frames
- **Transforms**: Define initial positioning and orientation
See also: [`Point`](@ref), [`Segment`](@ref), [`Group`](@ref), [`Pulley`](@ref)

# Arguments
- `name::String`: Physical model identifier (typically "ram" or "simple_ram")
- `set::Settings`: Configuration object containing simulation parameters

# Keyword Arguments
- `points::Vector{Point}=Point[]`: Point masses in the system
- `groups::Vector{Group}=Group[]`: Groups for wing deformation modeling
- `segments::Vector{Segment}=Segment[]`: Elastic connections between points
- `pulleys::Vector{Pulley}=Pulley[]`: Elements that redistribute line lengths
- `tethers::Vector{Tether}=Tether[]`: Groups of segments forming complete tether lines
- `winches::Vector{Winch}=Winch[]`: Ground-based winches for system control
- `wings::Vector{Wing}=Wing[]`: Rigid wing bodies
- `transforms::Vector{Transform}=Transform[]`: Initial positioning transformations

# Validation and Initialization
The constructor performs extensive validation:
- Ensures sequential indexing for all components
- Updates Settings object with component parameters
- Initializes geometric properties and unstretched lengths
- Applies transforms for initial positioning

# Returns
- `SystemStructure`: A complete system model ready for simulation

# Example
Create a minimal system with basic components:
  points = [Point(1, [0,0,0], STATIC), Point(2, [0,0,10], DYNAMIC)]
  segments = [Segment(1, (1,2), BRIDLE)]
  sys_struct = SystemStructure("test_system", set; points, segments)

Create a complete ram air wing system:
  wing = RamAirWing(set)
  sys_struct = SystemStructure(set, wing)  # Uses predefined ram model
"""
function SystemStructure(name, set; 
        points=Point[], 
        groups=Group[], 
        segments=Segment[], 
        pulleys=Pulley[], 
        tethers=Tether[], 
        winches=Winch[], 
        wings=Wing[],
        transforms=Transform[],
    )
    for (i, point) in enumerate(points)
        @assert point.idx == i
        @assert point.transform_idx <= length(transforms)
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
        set.l_tethers[i]   = winch.tether_length
        set.v_reel_outs[i] = winch.tether_vel
    end
    for (i, wing) in enumerate(wings)
        @assert wing.idx == i
    end
    for (i, transform) in enumerate(transforms)
        @assert transform.idx == i
        set.elevations[i] = rad2deg(transform.elevation)
        set.azimuths[i]   = rad2deg(transform.azimuth)
        set.headings[i]   = rad2deg(transform.heading)
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
    set.physical_model = name
    sys_struct = SystemStructure(name, set, points, groups, segments, pulleys, tethers, winches, wings, transforms, y, x, jac)
    init!(sys_struct, set)
    return sys_struct
end

function SystemStructure(set::Settings, wing::RamAirWing)
    length(set.bridle_fracs) != 4 && throw(ArgumentError("4 bridle fracs should be provided for all models."))

    if set.physical_model == "ram"
        return create_ram_sys_struct(set, wing)
    elseif set.physical_model == "simple_ram"
        return create_simple_ram_sys_struct(set, wing)
    else
        throw(ArgumentError("Undefined physical model"))
    end
end

function init!(transforms::Vector{Transform}, sys_struct::SystemStructure)
    @unpack points, wings = sys_struct
    for transform in transforms
        # ==================== TRANSLATE ==================== #
        base_pos, curr_base_pos = get_base_pos(transform, wings, points)
        T = base_pos - curr_base_pos
        for point in points
            if point.transform_idx == transform.idx
                point.pos_w .= point.pos_cad .+ T
            end
        end
        for wing in wings
            if wing.transform_idx == transform.idx
                wing.pos_w .= wing.pos_cad .+ T
            end
        end

        # ==================== ROTATE ==================== #
        curr_rot_pos = get_rot_pos(transform, wings, points)
        curr_elevation = KiteUtils.calc_elevation(curr_rot_pos - base_pos)
        curr_azimuth = -KiteUtils.azimuth_east(curr_rot_pos - base_pos)
        curr_R_t_w = calc_R_t_w(curr_elevation, curr_azimuth)
        R_t_w = calc_R_t_w(transform.elevation, transform.azimuth)

        for point in points
            if point.transform_idx == transform.idx
                vec = point.pos_w - base_pos
                vec_along_z = rotate_around_z(curr_R_t_w' * vec, transform.heading)
                point.pos_w .= base_pos + R_t_w * vec_along_z
            end
            if point.type == WING
                wing = wings[point.wing_idx]
                point.pos_b .= wing.R_b_c' * (point.pos_cad - wing.pos_cad) # TODO: test this
            end
        end
        for wing in wings
            if wing.transform_idx == transform.idx
                vec = wing.pos_w - base_pos
                vec_along_z = rotate_around_x(curr_R_t_w' * vec, transform.heading)
                wing.pos_w .= base_pos + R_t_w * vec_along_z
                for i in 1:3
                    wing.R_b_c[:, i] .= R_t_w * rotate_around_x(curr_R_t_w' * wing.R_b_c[:, i], transform.heading)
                end
            end
        end

    end
end

function calc_pos(wing::RamAirWing, gamma, frac)
    le_pos = [wing.le_interp[i](gamma) for i in 1:3]
    chord = [wing.te_interp[i](gamma) for i in 1:3] .- le_pos
    pos = le_pos .+ chord .* frac
    return pos
end

function create_tether(tether_idx, set, points, segments, tethers, attach_point, type, dynamics_type, z=[0,0,1])
    winch_pos = find_axis_point(attach_point.pos_cad, set.l_tether, z)
    dir = winch_pos - attach_point.pos_cad
    segment_idxs = Int16[]
    for i in 1:set.segments
        frac = i / set.segments
        pos = attach_point.pos_cad + frac * dir
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

function create_ram_sys_struct(set::Settings, vsm_wing::RamAirWing)
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
        points = [
            points
            Point(7+i_pnt, bridle_top[1], dynamics_type)
            Point(8+i_pnt, bridle_top[2], WING)
            Point(9+i_pnt, bridle_top[3], dynamics_type)
            Point(10+i_pnt, bridle_top[4], dynamics_type)

            Point(11+i_pnt, bridle_top[2] - 1z, dynamics_type)

            Point(12+i_pnt, bridle_top[1] - 2z, dynamics_type)
            Point(13+i_pnt, bridle_top[3] - 2z, dynamics_type)

            Point(14+i_pnt, bridle_top[1] - 4z, dynamics_type)
            Point(15+i_pnt, bridle_top[3] - 4z, dynamics_type)
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

    wings = [Wing(1, [1,2,3,4], I(3), zeros(3))]
    transforms = [Transform(1, deg2rad(set.elevation), deg2rad(set.azimuth), deg2rad(set.heading), zeros(3), points[end].idx; wing_idx=1)]
    
    return SystemStructure(set.physical_model, set; points, groups, segments, pulleys, tethers, winches, wings, transforms)
end

function create_simple_ram_sys_struct(set::Settings, wing::RamAirWing)
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

    wings = [Wing(1, [1,2,3,4], I(3), zeros(3))]
    transforms = [Transform(1, deg2rad(set.elevation), deg2rad(set.azimuth), deg2rad(set.heading), zeros(3), points[end].idx; wing_idx=1)]

    return SystemStructure(set.physical_model, set; points, groups, segments, pulleys, tethers, winches, wings, transforms)
end

function init!(sys_struct::SystemStructure, set::Settings)
    @unpack points, groups, segments, pulleys, tethers, winches, wings, transforms = sys_struct

    for segment in segments
        (segment.type === BRIDLE) && (segment.diameter = 0.001set.bridle_tether_diameter)
        (segment.type === POWER) && (segment.diameter = 0.001set.power_tether_diameter)
        (segment.type === STEERING) && (segment.diameter = 0.001set.steering_tether_diameter)
        @assert (0 < segment.diameter < 1)
    end

    for winch in winches
        winch.tether_length = set.l_tethers[winch.idx]
        winch.tether_vel    = set.v_reel_outs[winch.idx]
    end

    (length(groups) > 0) && (first_moment_frac = groups[1].moment_frac)
    for group in groups
        group.twist = 0.0
        group.twist_vel = 0.0
        @assert group.moment_frac ≈ first_moment_frac "All group.moment_frac must be the same."
    end
    
    for transform in transforms
        transform.elevation = deg2rad(set.elevations[transform.idx])
        transform.azimuth   = deg2rad(set.azimuths[transform.idx])
        transform.heading   = deg2rad(set.headings[transform.idx])
    end

    for segment in segments
        (segment.l0 ≈ 0) && (segment.l0 = norm(points[segment.point_idxs[1]].pos_cad - points[segment.point_idxs[2]].pos_cad))
        @assert (segment.l0 > 0)
    end

    for pulley in pulleys
        segment1, segment2 = segments[pulley.segment_idxs[1]], segments[pulley.segment_idxs[2]]
        pulley.sum_length = segment1.l0 + segment2.l0
        pulley.length = segment1.l0
        pulley.vel = 0.0
        @assert !(pulley.sum_length ≈ 0)
    end

    init!(transforms, sys_struct)
    return nothing
end
