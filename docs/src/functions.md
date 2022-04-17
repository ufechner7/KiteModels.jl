```@meta
CurrentModule = KiteModels
```
## Introduction
Most of the functions work on a KPS3 or KPS4 object. For this, the variable s is used.
Such a variable can be created with the lines:
```julia
using KiteModels, KitePodModels, KiteUtils
const s = KPS3(KCU(se()))
```
Or, if you want to use the 4 point kite model:
```julia
using KiteModels, KitePodModels, KiteUtils
const s = KPS4(KCU(se()))
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
calc_height
winch_force
spring_forces
lift_drag
lift_over_drag
v_wind_kite
```

## High level simulation interface
```@docs
init_sim
next_step
```

## Low level simulation interface
```@docs
clear
find_steady_state
residual!
```

## Environment
```@docs
calc_rho
calc_wind_factor
```

## Helper functions
```@docs
calc_drag
calc_set_cl_cd
calc_aero_forces
calc_particle_forces
initial_kite_ref_frame
inner_loop
loop
get_particles
```
