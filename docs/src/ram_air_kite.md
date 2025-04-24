```@meta
CurrentModule = KiteModels
```
## Introduction
The RamAirKite is based on ModelingToolkit, which allows to define the differential algebraic equations in symbolic form.

## Types
```@docs
Point
Pulley
Segment
Tether
Winch
KitePointGroup
```

## Functions
```@docs
create_sys!
diff_eqs!
reinit!
scalar_eqs!
linear_vsm_eqs!
force_eqs!
```