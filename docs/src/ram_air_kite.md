```@meta
CurrentModule = KiteModels
```
## Introduction
The [`SymbolicAWEModel`](@ref) is based on ModelingToolkit, which allows to define the differential algebraic equations in symbolic form. It does not use a KCU. Instead, the kite is controlled from the ground (e.g. a ship) using three winches and four tethers.

## Public constructors
```@docs
Point(::Any, ::Any, ::Any)
Group(::Any, ::Any, ::RamAirWing, ::Any, ::Any, ::Any)
Group(::Any, ::Any, ::Any, ::Any, ::Any, ::Any, ::Any)
Segment(::Any, ::Any, ::Any)
Pulley(::Any, ::Any, ::Any)
Winch(::Any, ::Any, ::Any, ::Any)
Wing(::Any, ::Any, ::Any, ::Any)
Transform(::Any, ::Any, ::Any, ::Any)
SystemStructure(::Any, ::Any)
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
Tether
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