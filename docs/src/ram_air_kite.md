```@meta
CurrentModule = KiteModels
```
## Introduction
The [`RamAirKite`](@ref) is based on ModelingToolkit, which allows to define the differential algebraic equations in symbolic form. It does not use a KCU. Instead, the kite is controlled from the ground (e.g. a ship) using three winches and four tethers.

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
KitePointGroup
PointMassSystem
```

## Private functions
```@docs
create_sys!
diff_eqs!
reinit!
scalar_eqs!
linear_vsm_eqs!
force_eqs!
```