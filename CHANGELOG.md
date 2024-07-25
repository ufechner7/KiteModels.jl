# Changelog

### KiteModels v0.6.0 - 2024-07-25
#### Changed
- use new version of Winchmodels.jl which provides an additional, torque controlled winch
- add many new winch parameters to settings.jl
- rename `v_ro` to `set_speed` in function step()

### KiteModels v0.5.16 - 2024-06-25
#### Changed
- bump KiteUtils to version 0.6.16
- bump ControlPlots to version 0.0.12

### KiteModels v0.5.15 - 2024-06-18
#### Changed
- bump KiteUtils to version 0.6.12
- drop support for Julia 1.9

### KiteModels v0.5.14 - 2024-05-05
#### Changed
- replace Plots with ControlPlots in the examples

### KiteModels v0.5.13 - 2024-04-14
#### Changed
- use `rel_compr_stiffness` and `rel_damping` from settings.yaml 

### KiteModels v0.5.12 - 2024-04-14
#### Changed
- update KiteUtils to v0.6.7
- update Documenter to v1.0

#### Added
- add type `KPS4_3L`, which is now only a copy of `KPS4`, but shall implement a kite with the steering
  lines going to the ground

### KiteModels v0.5.11 - 2024-04-04
#### Added
- document the support for the `DImplicitEuler` solver, which is not very accurate,
  but because it is well known it can serve as a reference
- support changing `max_order` for the `DFBDF`
- further reduced the memory usage

### KiteModels v0.5.10 - 2024-04-03
#### Added

- it is now possible (and suggested) to use the DAE solver DFBDF.

This requires adding the following line to the settings.yaml file: 

    solver: "DFBDF"

The new solver is much faster (4x average, 1.8x worst case), has a lot less memory allocations (~ 50%) and is also much more stable in highly dynamic situations.

### KiteModels v0.5.8 - 2024-04-01

#### Added
- new, non-allocating function `update_sys_state!(ss::SysState, s::AKM, zoom=1.0)`

### KiteModels v0.5.7 - 2024-04-01

#### Changed
- improved performance by 10% by implementing custom `norm()` function for 3D vectors

### KiteModels v0.5.6 - 2024-03-30

#### Fixed
- fix the method `clear!(s::KPS4)` which failed for models with less than 6 tether segments

Simulations should work fine now for one to about 28 tether segments (no hard upper limit, but things become slow and the visualization ugly if you have too many segments).
