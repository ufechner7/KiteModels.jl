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

## Public constructors
```@docs
SystemStructure(::Any, ::Any)
SymbolicAWEModel(::Settings, ::SystemStructure, ::Vector{<:BodyAerodynamics}, ::Vector{<:VortexStepMethod.Solver})
SymbolicAWEModel(::Settings)
Point(::Any, ::Any, ::Any)
Group(::Any, ::Any, ::RamAirWing, ::Any, ::Any, ::Any)
Group(::Any, ::Any, ::Any, ::Any, ::Any, ::Any, ::Any)
Segment(::Any, ::Any, ::Any)
Pulley(::Any, ::Any, ::Any)
Tether
Winch(::Any, ::Any, ::Any, ::Any)
Wing(::Any, ::Any, ::Any, ::Any)
Transform(::Any, ::Any, ::Any, ::Any)
```

## Private enumerations
```@docs
SegmentType
DynamicsType
```

## Private types
```@docs
Point
Pulley
Segment
Winch
Group
Wing
Transform
SystemStructure
```

## Private functions
```@docs
wing_eqs!
reinit!
scalar_eqs!
linear_vsm_eqs!
force_eqs!
```