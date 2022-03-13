```@meta
CurrentModule = KiteModels
```
## Introduction
Most of the functions work on a KPS3 object. For this, the variable s is used.
Such a variable can be created with the lines:
```julia
using KiteModels, KitePodSimulator
const s = KPS3(KCU())
```

## Input functions
```@docs
set_v_reel_out
set_depower_steering
set_v_wind_ground
```

## Output functions
```@docs
unstretched_length
tether_length
winch_force
spring_forces
lift_drag
lift_over_drag
v_wind_kite
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

## Helper functions
```@docs
clear
find_steady_state
calc_drag
calc_set_cl_cd
```
