<!--
SPDX-FileCopyrightText: 2025 Uwe Fechner

SPDX-License-Identifier: MIT
-->

# KiteModels

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ufechner7.github.io/KiteModels.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ufechner7.github.io/KiteModels.jl/dev)
[![CI](https://github.com/ufechner7/KiteModels.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/ufechner7/KiteModels.jl/actions/workflows/CI.yml)
[![Coverage](https://codecov.io/gh/ufechner7/KiteModels.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ufechner7/KiteModels.jl)
[![DOI](https://zenodo.org/badge/443855286.svg)](https://zenodo.org/doi/10.5281/zenodo.13310253)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

## Kite power system models, consisting of tether and kite
The models have the following subcomponents, implemented in separate packages:
- AtmosphericModel from [AtmosphericModels](https://github.com/aenarete/AtmosphericModels.jl)
- WinchModel from [WinchModels](https://github.com/aenarete/WinchModels.jl) 
- KitePodModel from  [KitePodModels](https://github.com/aenarete/KitePodModels.jl)
- The aerodynamic forces and moments of some of the models are calculated using the package [VortexStepMethod](https://github.com/Albatross-Kite-Transport/VortexStepMethod.jl)

This package is part of Julia Kite Power Tools, which consists of the following packages:
<p align="center"><img src="https://github.com/ufechner7/KiteModels.jl/blob/main/docs/src/kite_power_tools.png" width="500" /></p>

## News
#### Work in progress
- a new 5-point model based on ModelingToolkit (MTK) is in development;  
  this will allow to create linearized models around any operation point and to do analysis in the frequency domain.
#### April 2025
- a new model `SymbolicAWEModel` was contributed, based on the package [VortexStepMethod](https://github.com/Albatross-Kite-Transport/VortexStepMethod.jl)
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

## What to install
If you want to run simulations and see the results in 3D, please install the meta package  [KiteSimulators](https://github.com/aenarete/KiteSimulators.jl) . If you are not interested in 3D visualization or control you can just install this package.

## Installation
If possible, install [Julia 1.11](https://ufechner7.github.io/2024/08/09/installing-julia-with-juliaup.html), if you haven't already. Julia 1.10 is still supported, but the performance is worse. On Linux, make sure that Python3 and Matplotlib are installed:
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
You can also run the ram-air-kite example like this:
```julia
include("examples/ram_air_kite.jl")
```
This might take two minutes. To speed up the model initialization, you can create a system image:
```bash
cd bin
./create_sys_image
```
If you now launch Julia with `./bin/run_julia` and then run the above example again, it should be about three
times faster.

## Advanced installation
If you intend to modify or extend the code, it is suggested to fork the `KiteModels.jl` repository and to check out your fork:
```bash
git clone https://github.com/USERNAME/KiteModels.jl
```
where USERNAME is your github username.
Then compile a system image:
```bash
cd KiteModels.jl/bin
./create_sys_image
```
If you now launch julia with:
```bash
cd ..
./bin/run_julia
```
You can run the examples with:
```julia
menu()
```
You can also run the ram-air-kite example like this:
```julia
include("examples/ram_air_kite.jl")
```

## One point model
This model assumes the kite to be a point mass. This is sufficient to model the aerodynamic forces, but the dynamic concerning the turning action of the kite is not realistic.
When combined with a controller for the turn rate it can be used to simulate a pumping kite power system with medium accuracy.

## Four point model
This model assumes the kite to consist of four-point masses with aerodynamic forces acting on points B, C and D. It reacts much more realistically than the one-point model because it has rotational inertia in every axis.
<p align="center"><img src="https://github.com/ufechner7/KiteModels.jl/raw/main/docs/src/4-point-kite.png" width="200" /></p>

## Ram air kite model
This model represents the kite as a deforming rigid body, with orientation governed by quaternion dynamics. Aerodynamics are computed using the Vortex Step Method. The kite is controlled from the ground via four tethers.

## Tether
The tether is modeled as point masses, connected by spring-damper elements. Aerodynamic drag is modeled realistically. When reeling out or in the unstreched length of the spring-damper elements
is varied. This does not translate into physics directly, but it avoids adding point masses at run-time, which would be even worse because it would introduce discontinuities. When using
Dyneema or similar high-strength materials for the tether the resulting system is very stiff which is a challenge for the solver.

## Further reading
The models KPS3 and KPS4 are described in detail in [Dynamic Model of a Pumping Kite Power System](http://arxiv.org/abs/1406.6218).

## Replaying log files
If you want to replay old flight log files in 2D and 3D to understand and explain better how kite power systems work, please have a look at [KiteViewer](https://github.com/ufechner7/KiteViewer) . How new log files can be created and replayed is explained in the documentation of [KiteSimulators](https://github.com/aenarete/KiteSimulators.jl) .

## Licence
This project is licensed under the MIT and the MPL-2.0 License. The documentation is licensed under the CC-BY-4.0 License. Please see the below `Copyright notice` in association with the licenses that can be found in the file [LICENSE](LICENSE) in this folder.

## Copyright notice
Technische Universiteit Delft hereby disclaims all copyright interest in the package “KiteModels.jl” (models for airborne wind energy systems) written by the Author(s).

Prof.dr. H.G.C. (Henri) Werij, Dean of Aerospace Engineering, Technische Universiteit Delft.

See the copyright notices in the source files, and the list of authors in [AUTHORS.md](AUTHORS.md).

## Donations
If you like this software, please consider donating to [Flood in Kenya](https://www.gofundme.com/f/climate-refugees-in-kenya) .

## See also
- [Research Fechner](https://research.tudelft.nl/en/publications/?search=Fechner+wind&pageSize=50&ordering=rating&descending=true) for the scientic background of this code
- The meta-package  [KiteSimulators](https://github.com/aenarete/KiteSimulators.jl)
- the package [KiteUtils](https://github.com/ufechner7/KiteUtils.jl)
- the packages [WinchModels](https://github.com/aenarete/WinchModels.jl) and [KitePodModels](https://github.com/aenarete/KitePodModels.jl) and [AtmosphericModels](https://github.com/aenarete/AtmosphericModels.jl)
- the packages [WinchControllers](https://github.com/OpenSourceAWE/WinchControllers.jl), [KiteControllers](https://github.com/aenarete/KiteControllers.jl) and [KiteViewers](https://github.com/aenarete/KiteViewers.jl)
- the [VortexStepMethod](https://github.com/Albatross-Kite-Transport/VortexStepMethod.jl)

**Documentation** [Stable Version](https://ufechner7.github.io/KiteModels.jl/stable) --- [Development Version](https://ufechner7.github.io/KiteModels.jl/dev)
