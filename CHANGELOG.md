<!--
SPDX-FileCopyrightText: 2025 Uwe Fechner, Bart van de Lint
SPDX-License-Identifier: MIT
-->
### unreleased
- rename init_sim! to init!
- removed the parameter `upwind_dir!` from `init!`; use set.upwind_dir instead. Careful: This is in degrees.
- add the test script test_interface.jl
- add the field `integrator` to KPS4 and KPS3 structs

### KiteModels v0.8.1 2025-06-20
#### CHANGED
- renamed POWER to POWER_LINE and STEERING to STEERING_LINE
- improved documentation, fixed example

### KiteModels v0.8.0 2025-06-19
#### Added
- add a tutorial for custom system structures
- add documentation for `SystemStructure` components [#229](https://github.com/ufechner7/KiteModels.jl/pull/229)
- add the `Transform` struct which defines initial orientation of the `SystemStructure` [#214](https://github.com/ufechner7/KiteModels.jl/pull/214)
- add a custom amount of kites to the `SystemStructure` [#208](https://github.com/ufechner7/KiteModels.jl/pull/208)
- make model initialization faster [#222](https://github.com/ufechner7/KiteModels.jl/pull/222)
- implement missing methods for the `SymbolicAWEModel` [#198](https://github.com/ufechner7/KiteModels.jl/pull/198)
#### Changed
- fixed the performance regression when using the `KPS4` model with a winch model of type  `AsyncMachine`
- improved the script `create_sys_image`; it is now also available if you install the package without using git
- make model initialization faster [#222](https://github.com/ufechner7/KiteModels.jl/pull/222)
- make next_step! return nothing [#213](https://github.com/ufechner7/KiteModels.jl/pull/213)
- Change the names of `RamAirKite` to `SymbolicAWESystem` and `PointMassSystem` to `SystemStructure` [#208](https://github.com/ufechner7/KiteModels.jl/pull/208)

### KiteModels v0.7.4 2025-06-08
#### Added
- added licenses to each file, the command `pipx run reuse lint` succeeds now
- add the command above to the CI scripts
- the script `create_xz_file`
- the option to linearize the SymbolicAWEModel system using ModelingToolkit
- a simplified ram air kite model for faster development and testing
- the example `lin_ram_model.jl` to show how to linearize a model
- add the page `Examples SymbolicAWEModel` do the documentation
#### Changed
- the example `ram_air_kite.jl` can now be run like this `SIMPLE=true; include("examples/ram_air_kite.jl")`
- the package `Rotations` is no longer re-exported
- improved documentation
#### Fixed
- small fixes of the SymbolicAWEModel model

### KiteModels v0.7.3 2025-05-05
#### Fixed
- fix function update_sys_state!()

### KiteModels v0.7.2 2025-05-05
#### Changed
- bump KiteUtils to v0.10.5, which provides much more fields for in the SysState

### KiteModels v0.7.1 2025-04-28
#### Changed
- fixed or documented issues found by `Aqua.jl`
- made `DSP` a test dependency
- remove package `OrdinaryDiffEqSDIRK`
- improve documentation for `SymbolicAWEModel`
#### Added
- the examples `calc_spectrum.jl` and `plot_spectrum.jl` to the menu
- the quality insurance package `Aqua.jl`
- added the script `update_default_manifest`
- calculate `side_slip` angle in radian

### KiteModels v0.7.0 2025-04-20
#### Fixed
- fixed broken installation by freezing NLSolversBase to `~7.8.3` in Project.toml
#### Added
- added `mwe_26.jl` for debugging the initial state solver 
- the example `ram_air_kite.jl`
- the struct `SystemStructure` for easy definition of the kite power system
#### Changed
- BREAKING: the model KPS_3L was renamed to SymbolicAWEModel
- the SymbolicAWEModel model is using the **VortexStepMethod** with a deforming wing now
- bump KiteUtils to `v0.10`
- bump ModellingToolkit to `9.72`
- bump VortexStepMethod to `1.2.5`
- the file CONTRIBUTING.md was updated

### KiteModels v0.6.17 2025-02-11
#### Changed
- always use the brake if the `set_speed` is zero; this fixes the example `steering_test_4p.jl`

### KiteModels v0.6.16 2025-02-06
#### Changed
- `initial_reel_out_4p.jl` shows a simulation that starts with an initial reel-out speed > 0
- `initial_reel_out_4p_torque_control` runs a simulation with a torque controlled winch and an initial reel-out speed > 0
#### Fixed
- initial reel-out speed handled correctly

### KiteModels v0.6.15 2025-02-03
- log kcu_steering in SysState (output of KCU without applying corrections)
- fix tests for Julia 1.11.3
- cleanup `run_julia`
- add function `calculate_rotational_inertia!()`
- add example `calculate_rotational_inertia.jl` and add it to the menu

### KiteModels v0.6.14 2025-01-16
#### Fixed
- crash due to a new version of `DierckX_jll`
#### Changed
- bump versions in Project.toml
- add upwind_dir to `settings.yaml` file, remove `v_wind_ref` vector
- add `p_speed`, `i_speed` and `max_acc` to winch settings
- replaced `autodiff=false` with `autodiff=AutoFiniteDiff()` to fix warnings
- added package `ADTypes` to provide `AutoFiniteDiff()`
- cleanup code

### KiteModels v0.6.13 2024-12-06
#### Changed
- update the fields `set_steering`, `bearing` and `attractor` of the `SysState` struct 
  in the function `update_sys_state!`
- add the parameters `bearing` and `attractor` to the function `next_step!` for logging
#### Fixed
- fix #88: the function `init_sim!()` has the new parameter `upwind_dir` to define the
  initial wind direction 

### KiteModels v0.6.12 2024-12-01
#### Changed
- update the fields `set_torque`, `set_force`, `set_speed`, `alpha3`, `alpha4`, `roll`, `pitch`, `yaw`
  of the `SysState` struct in the function `update_sys_state!`
- add the parameter `set_force` to the function `next_step!` for logging
- the four point kite model KPS4 was extended to include aerodynamic damping of pitch oscillations;
  for this purpose, the parameters `cmq` and `cord_length` must be defined in `settings.yaml`
- the four point kite model KPS4 was extended to include the impact of the deformation of the
  kite on the turn rate; for this, the parameter `smc` must be defined in `settings.yaml`
- improve examples
- add the packages `JLD2` and `Colors` to the system image
#### Added
- add examples `calc_spectrum.jl` and `plot_spectrum` to plot the eigenfrequencies of the system
- the script `calculate_rotational_inertia.jl` for calculating the inertia matrix of the kite
- function `menu2()` which displays a menu with scripts for model verification

### KiteModels v0.6.11 2024-11-09
#### Fixed
- fixed bug in spring_forces(), it used 4000N hardcoded max force
#### Changed
- bump KiteUtils to 0.9.0
- the fields AoA, CL2, CD2 and the vectors v_wind_gnd, v_wind_200m and v_wind_kite are now updated 
  when converting KPS3 or KPS4 to SysState
#### Added
- add example test_steady_state.jl

### KiteModels v0.6.10 2024-11-01
- fixed the installation of the examples
- updated the documentation
- the package Rotations is now re-exported by KiteModels
- reduce number of dependencies of the examples

### KiteModels v0.6.9 2024-11-01
- added tests for calc_azimuth(s::AKM), the azimuth in wind reference frame
- re-enable logging of the angles of attack of the three plates
- `steering_test_4p.jl` now calculates both `c1` and `c2` of the turn-rate law
- the environment variable `NO_MTK` disables the pre-compilation of the `SymbolicAWEModel` model
  to save time during development
- the script `menu2.jl` for model verification was added

### KiteModels v0.6.8 - 2024-10-23
#### Changed
- the sign of the steering signal was changed. Now, a positive steering signal causes a positive turn rate.
  The turn rate is the derivative of the heading angle.
- all orientation tests pass now (calculation of roll, pitch, yaw, azimuth_north, elevation, heading)
- add example `steering_test_1p.jl`
- improve `steering_test_4p.jl`, use fully powered kite now
- the logged steering signal is now divided by `set.cs_4p`, because the new version of KitePodModels.jl multiplies the steering value with this constant
- update documentation regarding `steering` and `heading`

### KiteModels v0.6.7 - 2024-10-20
#### Changed
- renamed test_init.jl to test_init_4p.jl
- by default, `azimuth` in wind reference frame is now used
- the orientation is now represented in NED reference frame
- the method `calc_heading` has two new, optional parameters: neg_azimuth=false, one_point=false
- the definition of heading and azimuth has changed, which will require adaptions in the controller
#### Added
- example `plot_side_cl.jl`
- example `plot_cl_cd_plate.jl`
- example `steering_test.jl`
- example `test_init_1p.jl`
- the test script `test_orientation.jl` was contributed by Daan van Wolffelaar; currently, some of these tests are still broken (error of about 2%)
#### Fixed
- many of the examples; all examples of `menu.jl` now work

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
- the SymbolicAWEModel model was replaced by the pure ModelingToolkit (MTK) based version. This allows not only a much faster simulation, but the results are also much more accurate.

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
- add type `SymbolicAWEModel`, which is now only a copy of `KPS4`, but shall implement a kite with the steering
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
