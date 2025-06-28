
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
    SegmentType `POWER_LINE` `STEERING_LINE` `BRIDLE`

Type of segment.

# Elements
- POWER_LINE: Belongs to a power line
- STEERING_LINE: Belongs to a steering line
- BRIDLE: Belongs to the bridle
"""
@enum SegmentType begin
    POWER_LINE
    STEERING_LINE
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

A point mass.

$(TYPEDFIELDS)
"""
mutable struct Point
    const idx::Int16
    const transform_idx::Int16 # idx of wing used for initial orientation
    const wing_idx::Int16
    const pos_cad::KVec3
    const pos_b::KVec3 # pos relative to wing COM in body frame
    const pos_w::KVec3 # pos in world frame
    const vel_w::KVec3 # vel in world frame
    const type::DynamicsType
    mass::SimFloat
end

"""
    Point(idx, pos_cad, type; wing_idx=1, vel_w=zeros(KVec3), transform_idx=1, mass=0.0)

Constructs a Point object. A point can be of four different [`DynamicsType`](@ref)s:
- `STATIC`: the point doesn't move. ``\\ddot{\\mathbf{r}} = \\mathbf{0}``
- `DYNAMIC`: the point moves according to Newton's second law. ``\\ddot{\\mathbf{r}} = \\mathbf{F}/m``
- `QUASI_STATIC`: the acceleration is constrained to be zero, by solving a nonlinear problem. ``\\mathbf{F}/m = \\mathbf{0}``
- `WING`: the point has a static position in the rigid body wing frame. ``\\mathbf{r}_w = \\mathbf{r}_{wing} + \\mathbf{R}_{b\\rightarrow w} \\mathbf{r}_b``

where:
- ``\\mathbf{r}`` is the point position vector
- ``\\mathbf{F}`` is the net force acting on the point
- ``m`` is the point mass
- ``\\mathbf{r}_w`` is the position in world frame
- ``\\mathbf{r}_{wing}`` is the wing center position
- ``\\mathbf{R}_{b\\rightarrow w}`` is the rotation matrix from body to world frame
- ``\\mathbf{r}_b`` is the position in body frame

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
function Point(idx, pos_cad, type; wing_idx=1, vel_w=zeros(KVec3), transform_idx=1, mass=0.0)
    Point(idx, transform_idx, wing_idx, pos_cad, zeros(KVec3), zeros(KVec3), vel_w, type, mass)
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
All points within a group undergo the same twist rotation about the chord vector.

The governing equation is:
```math
\\begin{aligned}
\\tau = \\underbrace{\\sum_{i=1}^{4} r_{b,i} \\times (\\mathbf{F}_{b,i} \\cdot \\hat{\\mathbf{z}})}_{\\text{bridles}} + \\underbrace{r_a \\times (\\mathbf{F}_a \\cdot \\hat{\\mathbf{z}})}_{\\text{aero}}
\\end{aligned}
```

![System Overview](group_slice.svg)

where:
- ``\\tau`` is the total torque about the twist axis
- ``r_{b,i}`` is the position vector of bridle point ``i`` relative to the twist center
- ``\\mathbf{F}_{b,i}`` is the force at bridle point ``i``
- ``\\hat{\\mathbf{z}}`` is the unit vector along the twist axis (chord direction)
- ``r_a`` is the position vector of the aerodynamic center relative to the twist center
- ``\\mathbf{F}_a`` is the aerodynamic force at the group's aerodynamic center

The group can have two [`DynamicsType`](@ref)s:
- `DYNAMIC`: the group rotates according to Newton's second law: ``I\\ddot{\\theta} = \\tau``
- `QUASI_STATIC`: the rotational acceleration is zero: ``\\tau = 0``

# Arguments
- `idx::Int16`: Unique identifier for the group
- `point_idxs::Vector{Int16}`: Indices of points that move together with this group's twist
- `vsm_wing::RamAirWing`: Wing geometry object used to extract local chord and spanwise vectors
- `gamma`: Spanwise parameter (typically -1 to 1) defining the group's location along the wing
- `type::DynamicsType`: Dynamics type (DYNAMIC for time-varying twist, QUASI_STATIC for equilibrium)
- `moment_frac::SimFloat`: Chordwise position (0=leading edge, 1=trailing edge) about which the group rotates

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
a common twist deformation. See: [`Group(::Any, ::Any, ::RamAirWing, ::Any, ::Any, ::Any)`](@ref).

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

Constructs a Segment object representing an elastic spring-damper connection between two points.

The segment follows Hooke's law with damping and aerodynamic drag:

**Spring-Damper Force:**
```math
\\mathbf{F}_{spring} = \\left[k(l - l_0) - c\\dot{l}\\right]\\hat{\\mathbf{u}}
```

**Aerodynamic Drag:**
```math
\\mathbf{F}_{drag} = \\frac{1}{2}\\rho C_d A |\\mathbf{v}_a| \\mathbf{v}_{a,\\perp}
```

**Total Force:**
```math
\\mathbf{F}_{total} = \\mathbf{F}_{spring} + \\mathbf{F}_{drag}
```

where:
- ``k = \\frac{E \\pi d^2/4}{l}`` is the axial stiffness
- ``l`` is current length, ``l_0`` is unstretched length
- ``c = \\frac{\\xi}{c_{spring}} k`` is damping coefficient
- ``\\hat{\\mathbf{u}} = \\frac{\\mathbf{r}_2 - \\mathbf{r}_1}{l}`` is unit vector along segment
- ``\\dot{l} = (\\mathbf{v}_1 - \\mathbf{v}_2) \\cdot \\hat{\\mathbf{u}}`` is extension rate
- ``\\mathbf{v}_{a,\\perp}`` is apparent wind velocity perpendicular to segment

# Arguments
- `idx::Int16`: Unique identifier for the segment.
- `point_idxs::Tuple{Int16, Int16}`: Tuple containing the indices of the two points connected by this segment.
- `type::SegmentType`: Type of the segment (POWER_LINE, STEERING_LINE, BRIDLE).

# Keyword Arguments
- `l0::SimFloat=zero(SimFloat)`: Unstretched length of the segment. Calculated from point positions if zero.
- `compression_frac::SimFloat=0.1`: Compression fraction of stiffness for compression behavior.

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

Constructs a Pulley object that enforces length redistribution between two segments.

The pulley constraint maintains constant total length while allowing force transmission:

**Constraint Equations:**
```math
l_1 + l_2 = l_{total} = \\text{constant}
```

**Force Balance:**
```math
F_{pulley} = F_1 - F_2
```

**Dynamics:**
```math
m\\ddot{l}_1 = F_{pulley} = F_1 - F_2
```

where:
- ``l_1, l_2`` are the lengths of connected segments
- ``F_1, F_2`` are the spring forces in the segments  
- ``m = \\rho_{tether} \\pi (d/2)^2 l_{total}`` is the total mass of both segments
- ``\\dot{l}_1 + \\dot{l}_2 = 0`` (velocity constraint)

The pulley can have two [`DynamicsType`](@ref)s:
- `DYNAMIC`: the length redistribution follows Newton's second law: ``m\\ddot{l}_1 = F_1 - F_2``
- `QUASI_STATIC`: the forces are balanced instantaneously: ``F_1 = F_2``

# Arguments
- `idx::Int16`: Unique identifier for the pulley.
- `segment_idxs::Tuple{Int16, Int16}`: Tuple containing the indices of the two segments connected by this pulley.
- `type::DynamicsType`: Dynamics type of the pulley (DYNAMIC or QUASI_STATIC).

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
    Tether(idx, segment_idxs)

Constructs a Tether object representing a flexible line composed of multiple segments.

A tether enforces a shared unstretched length constraint across all its constituent segments:

**Length Constraint:**
```math
\\sum_{i \\in \\text{segments}} l_{0,i} = L
```

**Winch Control:**
The unstretched tether length is controlled by winch acceleration:
```math
\\ddot L = \\alpha(v, F, u)
```

where:
- ``L`` is the tether length
- ``l_{0,i}`` is the segment unstretched length
- ``\\alpha(v, F, u)`` is the winch acceleration function depending on model type

# Arguments
- `idx::Int16`: Unique identifier for the tether
- `segment_idxs::Vector{Int16}`: Indices of segments that form this tether

# Returns
- `Tether`: A new Tether object

# Example
Create a tether from segments 1, 2, and 3:
```julia
    tether = Tether(1, [1, 2, 3])
```
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
    tether_length::Union{SimFloat, Nothing}
    tether_vel::SimFloat
end

"""
    Winch(idx, model, tether_idxs; tether_length=nothing, tether_vel=0.0)

Constructs a Winch object that controls tether length through torque or speed regulation.

**Tether Length Control:**
```math
\\ddot{L} = \\alpha(v, F, u)
```
where:
- ``L`` is the tether length
- ``v`` is the reel out velocity (tether extension rate)
- ``F`` is the tether force
- ``u`` is the applied torque or speed setpoint
- ``\\alpha(v, F, u)`` is the winch acceleration function depending on model type

where the winch acceleration function `f_winch` depends on the winch model type:
- **Torque-controlled**: Direct torque input with motor dynamics
- **Speed-controlled**: Velocity regulation with internal control loops

For detailed mathematical models of winch dynamics, motor characteristics, and control algorithms, 
see the [WinchModels.jl documentation](https://github.com/aenarete/WinchModels.jl/blob/main/docs/winch.md).

# Arguments
- `idx::Int16`: Unique identifier for the winch.
- `model::AbstractWinchModel`: The winch model (TorqueControlledMachine, AsyncMachine, etc.).
- `tether_idxs::Vector{Int16}`: Vector containing the indices of the tethers connected to this winch.

# Keyword Arguments
- `tether_vel::SimFloat=0.0`: Initial tether velocity (reel-out rate).
- `tether_length::SimFloat`: Initial tether length.

# Returns
- `Winch`: A new Winch object.

# Example
To create a Winch:
```julia
    winch = Winch(1, TorqueControlledMachine(set), [1, 2], 100.0)
```
"""
function Winch(idx, model, tether_idxs; tether_length=0.0, tether_vel=0.0)
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
    R_b_w::Matrix{SimFloat}
    angular_vel::KVec3
    pos_w::KVec3
    pos_cad::KVec3
    vel_w::KVec3
end
function Base.getproperty(wing::Wing, sym::Symbol)
    if sym == :orient
        return rotation_matrix_to_quaternion(wing.R_b_w)
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

**Rigid Body Dynamics:**
The wing follows standard rigid body equations of motion:

```math
\\begin{aligned}
\\frac{\\delta \\mathbf{q}_b^w}{\\delta t} &= \\frac{1}{2} \\Omega(\\boldsymbol{\\omega}_b) \\mathbf{q}_b^w \\\\
\\boldsymbol{\\tau}_w &= \\mathbf{I} \\frac{\\delta \\boldsymbol{\\omega}}{\\delta t} + \\boldsymbol{\\omega}_b \\times (\\mathbf{I}\\boldsymbol{\\omega}_b)
\\end{aligned}
```

where:
- ``\\mathbf{q}_b^w`` is the quaternion from world to body frame
- ``\\boldsymbol{\\omega}_b`` is the angular velocity in body frame
- ``\\Omega(\\boldsymbol{\\omega}_b)`` is the quaternion multiplication matrix
- ``\\mathbf{I}`` is the inertia tensor in body frame
- ``\\boldsymbol{\\tau}_w`` is the total applied torque to the rigid wing body (aerodynamic + tether forces)

**Coordinate Transformations:**
Points attached to the wing transform as:
```math
\\mathbf{r}_w = \\mathbf{r}_{w} + \\mathbf{R}_{b \\rightarrow w} \\mathbf{r}_b
```

where:
- ``\\mathbf{r}_w`` is the position in world frame
- ``\\mathbf{r}_{w}`` is the wing position in world frame
- ``\\mathbf{R}_{b \\rightarrow w}`` is the rotation from body to world frame
- ``\\mathbf{r}_b`` is the point position in body frame

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
    return Wing(idx, group_idxs, transform_idx, R_b_c, zeros(3,3), angular_vel, pos_w, pos_cad, vel_w)
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

Constructs a Transform object that orients system components using spherical coordinates.

**All points and wings with matching `transform_idx` are transformed together as a rigid body:**
1. **Translation**: Translate such that base is at specified base pos
1. **Rotation 1**: Rotate so target is at (elevation, azimuth) relative to base
2. **Rotation 2**: Rotate all components by `heading` around the base-target vector

```math
\\mathbf{r}_{transformed} = \\mathbf{r}_{base} + \\mathbf{R}_{heading} \\circ \\mathbf{R}_{elevation,azimuth}(\\mathbf{r} - \\mathbf{r}_{base})
```

# Arguments
- `idx::Int16`: Unique identifier for the transform
- `elevation::SimFloat`: Target elevation angle from base (radians)
- `azimuth::SimFloat`: Target azimuth angle from base (radians)  
- `heading::SimFloat`: Rotation around base-target vector (radians)

# Keyword Arguments
**Base Reference (choose one):**
- `base_pos + base_point_idx`: Fixed position and reference point
- `base_transform_idx`: Chain to another transform's position

**Target Object (choose one):**
- `wing_idx`: Wing to position at (elevation, azimuth)
- `rot_point_idx`: Point to position at (elevation, azimuth)

# Returns
- `Transform`: Transform affecting all components with matching `transform_idx`

# Examples
```julia
# Position wing and all associated points at 45° elevation, 30° azimuth
transform = Transform(1, deg2rad(45), deg2rad(30), 0.0; 
                     base_pos=[0,0,0], base_point_idx=1, wing_idx=1)

# Chain transforms for multi-kite systems
transform2 = Transform(2, deg2rad(30), deg2rad(45), deg2rad(10); 
                      base_transform_idx=1, wing_idx=2)
```
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

A discrete mass-spring-damper representation of a kite system, where point masses 
connected by elastic segments model the kite and tether dynamics.

# Components
- [`Point`](@ref): Point masses representing wing attachment points, dynamic bridle/tether points, and fixed ground anchors
- [`Group`](@ref): Collections of points that move together according to wing deformation (twist and trailing edge deflection)
- [`Segment`](@ref): Spring-damper elements connecting points
- [`Pulley`](@ref): Elements that redistribute line lengths between segments
- [`Tether`](@ref): Groups of segments with a common unstretched length
- [`Winch`](@ref): Ground-based winches that control tether lengths
- [`Wing`](@ref): Rigid wing bodies that serve as reference frames
- [`Transform`](@ref): Spatial transformations for initial positioning and orientation

See the individual component documentation for detailed mathematical models and governing equations.
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

## Components

- **Points** - See [`Point`](@ref) for discrete mass dynamics
- **Segments** - See [`Segment`](@ref) for elastic spring-damper connections  
- **Groups** - See [`Group`](@ref) for wing twist deformation modeling
- **Wings** - See [`Wing`](@ref) for rigid body dynamics
- **Pulleys** - See [`Pulley`](@ref) for length redistribution between segments
- **Tethers** - See [`Tether`](@ref) for segment groups with shared unstretched length
- **Winches** - See [`Winch`](@ref) for ground-based tether length control
- **Transforms** - See [`Transform`](@ref) for initial positioning and orientation

## Physical Models
- **"ram"**: 4 deformable wing groups, complex pulley bridle system
- **"simple_ram"**: 4 deformable wing groups, direct bridle connections

# Arguments
- `name::String`: Model identifier. "ram" and "simple_ram" are defined inside KiteModels.jl, provide a different name for a custom model.
- `set::Settings`: Configuration parameters (see [KiteUtils.Settings](https://ufechner7.github.io/KiteUtils.jl/stable/types/#KiteUtils.Settings))

# Returns
- `SystemStructure`: Complete system ready for building a [`SymbolicAWEModel`](@ref)

# Examples
```julia
# Auto-generate from wing geometry
wing = RamAirWing(set)
sys_struct = SystemStructure(set, wing)

# Manual construction
points = [Point(1, [0,0,0], STATIC), Point(2, [0,0,10], DYNAMIC)]
segments = [Segment(1, (1,2), BRIDLE)]
sys_struct = SystemStructure("custom", set; points, segments)
```
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
        (segment.l0 ≈ 0) && (segment.l0 = norm(points[segment.point_idxs[1]].pos_cad - points[segment.point_idxs[2]].pos_cad))
    end
    for (i, pulley) in enumerate(pulleys)
        @assert pulley.idx == i
    end
    for (i, tether) in enumerate(tethers)
        @assert tether.idx == i
    end
    for (i, winch) in enumerate(winches)
        @assert winch.idx == i
        if iszero(winch.tether_length)
            for segment_idx in tethers[winch.tether_idxs[1]].segment_idxs
                winch.tether_length += segments[segment_idx].l0
            end
        end
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

function apply_heading(vec, R_t_w, curr_R_t_w, heading)
    vec_along_z = rotate_around_z(curr_R_t_w' * vec, heading)
    return R_t_w * vec_along_z
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
                point.pos_w .= base_pos + apply_heading(vec, R_t_w, curr_R_t_w, transform.heading)
            end
            if point.type == WING
                wing = wings[point.wing_idx]
                point.pos_b .= wing.R_b_c' * (point.pos_cad - wing.pos_cad) # TODO: test this
            end
        end
        for wing in wings
            if wing.transform_idx == transform.idx
                vec = wing.pos_w - base_pos
                wing.pos_w .= base_pos + apply_heading(vec, R_t_w, curr_R_t_w, transform.heading)
                for i in 1:3
                    wing.R_b_w[:, i] .= apply_heading(wing.R_b_c[:, i], R_t_w, curr_R_t_w, transform.heading)
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

    points, segments, tethers, left_power_idx = create_tether(1, set, points, segments, tethers, attach_points[1], POWER_LINE, dynamics_type, z)
    points, segments, tethers, right_power_idx = create_tether(2, set, points, segments, tethers, attach_points[3], POWER_LINE, dynamics_type, z)
    points, segments, tethers, left_steering_idx = create_tether(3, set, points, segments, tethers, attach_points[2], STEERING_LINE, dynamics_type, z)
    points, segments, tethers, right_steering_idx = create_tether(4, set, points, segments, tethers, attach_points[4], STEERING_LINE, dynamics_type, z)

    winches = [winches; Winch(1, TorqueControlledMachine(set), [left_power_idx, right_power_idx])]
    winches = [winches; Winch(2, TorqueControlledMachine(set), [left_steering_idx])]
    winches = [winches; Winch(3, TorqueControlledMachine(set), [right_steering_idx])]

    wings = [Wing(1, [1,2,3,4], I(3), zeros(3))]
    transforms = [Transform(1, deg2rad(set.elevation), deg2rad(set.azimuth), deg2rad(set.heading);
                                    base_pos= zeros(3), base_point_idx=points[end].idx, wing_idx=1)]
    
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

    points, segments, tethers, left_power_idx = create_tether(1, set, points, segments, tethers, points[5], POWER_LINE, dynamics_type)
    points, segments, tethers, right_power_idx = create_tether(2, set, points, segments, tethers, points[7], POWER_LINE, dynamics_type)
    points, segments, tethers, left_steering_idx = create_tether(3, set, points, segments, tethers, points[6], STEERING_LINE, dynamics_type)
    points, segments, tethers, right_steering_idx = create_tether(4, set, points, segments, tethers, points[8], STEERING_LINE, dynamics_type)

    winches = [winches; Winch(1, TorqueControlledMachine(set), [left_power_idx, right_power_idx], set.l_tether)]
    winches = [winches; Winch(2, TorqueControlledMachine(set), [left_steering_idx], set.l_tether)]
    winches = [winches; Winch(3, TorqueControlledMachine(set), [right_steering_idx], set.l_tether)]

    wings = [Wing(1, [1,2,3,4], I(3), zeros(3))]
    transforms = [Transform(1, deg2rad(set.elevation), deg2rad(set.azimuth), deg2rad(set.heading);
                                    base_pos= zeros(3), base_point_idx=points[end].idx, wing_idx=1)]

    return SystemStructure(set.physical_model, set; points, groups, segments, pulleys, tethers, winches, wings, transforms)
end

function init!(sys_struct::SystemStructure, set::Settings)
    @unpack points, groups, segments, pulleys, tethers, winches, wings, transforms = sys_struct

    for segment in segments
        (segment.type === BRIDLE) && (segment.diameter = 0.001set.bridle_tether_diameter)
        (segment.type === POWER_LINE) && (segment.diameter = 0.001set.power_tether_diameter)
        (segment.type === STEERING_LINE) && (segment.diameter = 0.001set.steering_tether_diameter)
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
