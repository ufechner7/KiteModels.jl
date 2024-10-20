# Changelog
### KiteModels v0.6.7 - 2024-10-20
#### Changed
- renamed test_init.jl to test_init_4p.jl
- by default, `azimuth` in wind reference frame is now used
- the method `calc_heading` has two new, optional parameters: neg_azimuth=false, one_point=false
- the definition of heading and azimuth has changed, which will require adaptions in the controller
#### Added
- example `plot_side_cl.jl`
- example `plot_cl_cd_plate.jl`
- example `steering_test.jl`
- example `test_init_1p.jl`
- the test script `test_orientation.jl`; currently, some of these tests still fail (error of about 2%)
#### Fixed
- many of the examples
- the orientation is now represented in NED reference frame

### KiteModels v0.6.6 - 2024-09-03
#### Changed
- the method `next_step!` uses now `upwind_dir` as parameter and not longer `wind_dir`
- install `matplotlib` if it is not already installed after user confirmation in a Julia specific environment
- replaced OrdinaryDiffEq with the three packages OrdinaryDiffEqCore, OrdinaryDiffEqBDF
  and OrdinaryDiffEqSDIRK. This should help to reduce the pre-compilation time.
- set the parameter delta in the examples
- always specify the `system.yaml` file to use in the examples, always use `load_settings` instead of `se`. 
This ensures that the settings are always freshly loaded from the file when the script is launched, so any changes 
to the settings become immediately effective.
- the KPS4_3L model was replaced by the pure ModelingToolkit (MTK) based version. This allows not only a much faster simulation, but the results are also much more accurate.

### KiteModels v0.6.5 - 2024-08-12
#### Changed
- bump KiteUtils to 0.7.7
- add new examples to menu
- major change to the function that finds the initial equilibrium; the function `init_sim!` has the new
  parameter `delta` which should be in the range of 0.01 to 0.03.
- better error message if `init_sim!`, but no exception any more. I just returns `nothing`.
- remove dependency StatProfilerHTML
#### Added
- add KCU drag, based on kcu_diameter and cd_kcu
- add function bridle_length (not exported)
- unit tests for the KPS3_3L model, based on ModelingToolkit
- script `examples/plot_cl_cd.jl`
- script `examples/plot_cl_cd_plate.jl`
- script `torque_controlled_mtk.jl`
#### Fixed
- correct tether drag based on l_bridle; if the kite has more than 7 bridle lines l_bridle must be larger than bridle_length(se)

### KiteModels v0.6.4 - 2024-08-12
#### Added
- a new kite model, KPS3_3L was contributed by Bart van de Lint. It uses three lines to the ground and three winches for steering a ram-air foil kite.
- caching for the initial equilibrium
- function `KiteModels.install_examples()` that copies the examples, the data folder and installs the required extra packages
- log alpha2, alpha3, alpha4; they must never become negative
#### Fixed
- the calculation of the call-backs per time step was fixed in all examples and the tests

### KiteModels v0.6.3 - 2024-08-06
#### Changed
- the function `copy_examples()` copies now all examples
- updated the documentation
- improved the scripts in the bin folder, not relevant for most users

### KiteModels v0.6.2 - 2024-08-06
#### Changed
- renamed the example `simulate.jl` to `simulate_simple.jl`
- renamed the example `simulate_ii.jl` to `simulate_steering.jl`
- add the script `menu.jl` that provides a menu with all the examples
- bump ControlPlots to 0.1.4
- bump KiteUtils to 0.7.4

### KiteModels v0.6.1 - 2024-07-25
#### Changed
- bump WinchModels to 0.3.2
- bump KitePodModels to 0.3.3
- fix example `reel_out_4p_torque_control.jl`

### KiteModels v0.6.0 - 2024-07-25
#### Changed
- use a new version of `WinchModels.jl` which provides an additional, torque-controlled winch
- add many new winch parameters to `settings.jl`
- BREAKING change: rename `v_ro` to `set_speed` in function step()

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
