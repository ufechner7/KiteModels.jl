```@meta
CurrentModule = KiteModels
```

# KiteModels
Documentation for the package [KiteModels](https://github.com/ufechner7/KiteModels.jl).

The models have the following subcomponents, implemented in separate packages:
- AtmosphericModel from [AtmosphericModels](https://github.com/aenarete/AtmosphericModels.jl)
- WinchModel from [WinchModels](https://github.com/aenarete/WinchModels.jl) 
- KitePodModel from  [KitePodModels](https://github.com/aenarete/KitePodModels.jl)
- The aerodynamic forces and moments of some of the models are calculated using the package [VortexStepMethod](https://github.com/Albatross-Kite-Transport/VortexStepMethod.jl)

This package is part of Julia Kite Power Tools, which consist of the following packages:

![Julia Kite Power Tools](kite_power_tools.png)

## What to install
If you want to run simulations and see the results in 3D, please install the meta package  [KiteSimulators](https://github.com/aenarete/KiteSimulators.jl) which contains all other packages. If you are not interested in 3D visualization or control you can just install this package. When you have installed the package KiteSimulators, use the command `using KiteSimulators` instead of `using KiteModels` when this is mentioned in the documentation.

## Installation
Install [Julia 1.10](https://ufechner7.github.io/2024/08/09/installing-julia-with-juliaup.html) or later, if you haven't already. On Linux, make sure that Python3 and Matplotlib are installed:
```
sudo apt install python3-matplotlib
```
Before installing this software it is suggested to create a new project, for example like this:
```bash
mkdir test
cd test
julia --project="."
```
Then add KiteModels from  Julia's package manager, by typing:
```julia
using Pkg
pkg"add KiteModels"
``` 
at the Julia prompt. You can run the unit tests with the command (careful, can take 60 min):
```julia
pkg"test KiteModels"
```
You can copy the examples to your project with:
```julia
using KiteModels
KiteModels.install_examples()
```
This also adds the extra packages, needed for the examples to the project. Furthermore, it creates a folder `data`
with some example input files. You can now run the examples with the command:
```julia
include("examples/menu.jl")
```

## News
#### Work in progress
- a new 5-point model based on ModelingToolkit (MTK) is in development;  
  this will allow to create linearized models around any operation point and to do analysis in the frequency domain.
#### April 2025
- a new model `RamAirKite` was contributed, based on the package [VortexStepMethod](https://github.com/Albatross-Kite-Transport/VortexStepMethod.jl)
#### November 2024
- the four point kite model KPS4 was extended to include aerodynamic damping of pitch oscillations;
  for this purpose, the parameters `cmq` and `cord_length` must be defined in `settings.yaml`
- the four point kite model KPS4 was extended to include the impact of the deformation of the
  kite on the turn rate; for this, the parameter `smc` must be defined in `settings.yaml`
#### October 2024
- the orientation is now represented with respect to the NED reference frame
- azimuth is now calculated in wind reference frame. This allows it to handle changes of the wind direction
  during flight correctly.
- many unit tests added by a new contributor
- many tests for model verification added; they can be accessed using the `menu2.jl` script
- the documentation was improved
#### August 2024
- a new kite model, KPS3_3L was contributed. It uses three lines to the ground and three winches for steering a ram-air foil kite.
- a first [ModelingToolkit](https://docs.sciml.ai/ModelingToolkit/stable/) based model was added, which shows a much better performance and easier to read code
- a new KCU model was added which assumes a linear relationship between the depower settings and the depower angle and thus is easier to configure than the original model.
- the drag of the KCU is now taken into account
- the drag of the bridle is now taken into account correctly, also if the real kite has more bridle lines than the model
- the function to find the initial state is now much more robust
#### July 2024
- a new groundstation / winch type is now supported, the `TorqueControlledMachine`. It can be configured in the section `winch` of the `settings.yaml` file. It uses a set torque as input.
- a Python interface is now provided, see: [pykitemodels](https://github.com/ufechner7/pykitemodels)
#### April 2024
- added support for the native Julia DAE solver DFBDF. It is much more accurate and faster than the IDA solver that was used before.

## Provides
The type [`AbstractKiteModel`](@ref) with the implementation [`KPS3`](@ref), [`KPS4`](@ref) and [`RamAirKite`](@ref), representing the one point, the four point kite model and the ram air kite model, together with the high level simulation interface consisting of the functions [`init_sim!`](@ref) and [`next_step!`](@ref). Other kite models can be added inside or outside of this package by implementing the non-generic methods required for an AbstractKiteModel.

Additional functions to provide inputs and outputs of the model on each time step. In particular the constructor [`SysState`](@ref) can be called once per time step to create a SysState struct for
logging or for displaying the state in a viewer. For the KPS3 and KPS4 model, once per time step the [`residual!`](@ref) function is called as many times as needed to find the solution at the end
of the time step. The formulas are based on basic physics and aerodynamics and can be quite simple because a differential algebraic notation is used.

## One point model
This model assumes the kite to be a point mass. This is sufficient to model the aerodynamic forces, but the dynamic concerning the turning action of the kite is not realistic.
When combined with a controller for the turn rate it can be used to simulate a pumping kite power system with medium accuracy.

## Four point model
This model assumes the kite to consist of four-point masses with aerodynamic forces acting on points B, C and D. It reacts much more realistically than the one-point model because it has rotational inertia in every axis.

![Four point kite power system model](kps4.png)

## Ram air kite model
This model represents the kite as a deforming rigid body, with orientation governed by quaternion dynamics. Aerodynamics are computed using the Vortex Step Method. The kite is controlled from the ground via four tethers.

## Tether
The tether is modeled as point masses, connected by spring-damper elements. Aerodynamic drag is modeled realistically. When reeling out or in the unstreched length of the spring-damper elements
is varied. This does not translate into physics directly, but it avoids adding point masses at run-time, which would be even worse because it would introduce discontinuities. When using
Dyneema or similar high-strength materials for the tether the resulting system is very stiff which is a challenge for the solver.

## Reference frames and control inputs
- a positive `set_torque` will accelerate the reel-out, a negative `set_torque` counteract the pulling force of the kite. The unit is [N/m] as seen at the motor/generator axis.
- the `depower` settings are dimensionless and can be between zero and one. A value equal to $\mathrm{depower\_zero}/100$ from the `settings.yaml` file means that the kite is fully powered. 
- the `heading` angle, the direction the nose of the kite is pointing to is positive in clockwise direction when seen from above.
- the `steering` input, dimensionless and in the range of -1.0 .. 1.0. A positive steering input causes a positive turn rate (derivative of the heading).

A definition of the reference frames can be found [here](https://ufechner7.github.io/KiteUtils.jl/dev/reference_frames/) .

## Further reading
The one point and four point kite models are described in detail in [Dynamic Model of a Pumping Kite Power System](http://arxiv.org/abs/1406.6218).

## See also
- [Research Fechner](https://research.tudelft.nl/en/publications/?search=Fechner+wind&pageSize=50&ordering=rating&descending=true) for the scientic background of this code
- The meta-package  [KiteSimulators](https://github.com/aenarete/KiteSimulators.jl)
- the package [KiteUtils](https://github.com/ufechner7/KiteUtils.jl)
- the packages [WinchModels](https://github.com/aenarete/WinchModels.jl) and [KitePodModels](https://github.com/aenarete/KitePodModels.jl) and [AtmosphericModels](https://github.com/aenarete/AtmosphericModels.jl)
- the packages [KiteControllers](https://github.com/aenarete/KiteControllers.jl) and [KiteViewers](https://github.com/aenarete/KiteViewers.jl)
- the [VortexStepMethod](https://github.com/Albatross-Kite-Transport/VortexStepMethod.jl)

Authors: Uwe Fechner (uwe.fechner.msc@gmail.com), Bart van de Lint (bart@vandelint.net)
