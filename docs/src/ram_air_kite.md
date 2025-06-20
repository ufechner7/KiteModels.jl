```@meta
CurrentModule = KiteModels
```
## Introduction
The [`SystemStructure`](@ref) provides a flexible framework for defining the physical structure of airborne wind energy (AWE) systems using discrete mass-spring-damper models. This structure can represent many different AWE system configurations, from simple single-line kites to complex multi-wing systems with intricate bridle networks.

The [`SystemStructure`](@ref) serves as input to the [`SymbolicAWEModel`](@ref), which is based on ModelingToolkit and automatically generates symbolic differential algebraic equations from the structural definition.

## Workflow
1. Define system components ([`Point`](@ref), [`Segment`](@ref), [`Group`](@ref), etc.) 
2. Assemble into a [`SystemStructure`](@ref)
3. Pass to [`SymbolicAWEModel`](@ref) for automatic MTK model generation
4. Simulate the resulting symbolic model

## Public enumerations
```@docs
SegmentType
DynamicsType
```

## Public constructors
```@docs
SystemStructure(name, set; points=Point[], groups=Group[], segments=Segment[], 
                   pulleys=Pulley[], tethers=Tether[], winches=Winch[], 
                   wings=Wing[], transforms=Transform[])
SystemStructure
SymbolicAWEModel(::Settings, ::SystemStructure, ::Vector{<:BodyAerodynamics}, ::Vector{<:VortexStepMethod.Solver})
SymbolicAWEModel(::Settings)
Point(idx, pos_cad, type)
Point
Group(::Any, ::Any, ::RamAirWing, ::Any, ::Any, ::Any)
Group
Segment(idx, point_idxs, type)
Segment
Pulley(idx, segment_idxs, type)
Pulley
Tether
Winch(idx, model, tether_idxs, tether_length; tether_vel=0.0)
Winch
Wing(idx, group_idxs, R_b_c, pos_cad)
Wing
Transform(idx, elevation, azimuth, heading)
Transform
```

## Private functions
```@docs
wing_eqs!
reinit!
scalar_eqs!
linear_vsm_eqs!
force_eqs!
```