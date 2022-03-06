```@meta
CurrentModule = KiteModels
```
## Introduction
Most of the functions work on a KPS3 object. For this, the variable s is used.
Such a variable can be created with the lines:
```julia
using KiteUtils, KitePodSimulator
const s = KPS3(KCU())
```

## Input functions
```@docs
set_v_reel_out
set_l_tether
set_depower_steering
```

## Output functions
```@docs
get_l_tether
get_force
get_spring_forces
get_lift_drag
get_lod
```

## Callback function for the DAE solver
```@docs
residual!
```

## Environment
```@docs
calc_rho
calc_wind_factor
```
