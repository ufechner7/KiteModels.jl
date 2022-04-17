```@meta
CurrentModule = KiteModels
```

# KiteModels

Documentation for the package [KiteModels](https://github.com/ufechner7/KiteModels.jl).

## Installation

Download [Julia 1.6](http://www.julialang.org) or later, if you haven't already. You can add KiteModels from  Julia's package manager, by typing 
```
] add KiteModels
``` 
at the Julia prompt.

If you are using Windows, it is suggested to install git and bash, too. This is explained for example here: [Julia on Windows](https://github.com/ufechner7/KiteViewer/blob/main/doc/Windows.md) .

## Provides

The type [`AbstractKiteModel`](@ref) with the implementation [`KPS3`](@ref) and [`KPS4`](@ref), representing the model, together with the high level simulation interface consisting of the functions
[`init_sim!`](@ref) and [`next_step!`](@ref). Other kite models can be added inside or outside of this package by implementing the non-generic methods required for an AbstractKiteModel.

Additional functions to provide inputs and outputs of the model on each time step. Per time step the [`residual!`](@ref) function is called as many times as needed to find the solution at the end
of the time step. The formulas are based on basic physics and aerodynamics and can be quite simple because a differential algebraic notation is used.

![Four point kite power system model](kps4.png)

## Further reading
These models are described in detail in [Dynamic Model of a Pumping Kite Power System](http://arxiv.org/abs/1406.6218).

## See also
- [Research Fechner](https://research.tudelft.nl/en/publications/?search=Uwe+Fechner&pageSize=50&ordering=rating&descending=true) for the scientic background of this code
- The application [KiteViewer](https://github.com/ufechner7/KiteViewer)
- the package [KiteUtils](https://github.com/ufechner7/KiteUtils.jl)
- the package [KitePodSimulator](https://github.com/ufechner7/KitePodSimulator.jl)

Author: Uwe Fechner (uwe.fechner.msc@gmail.com)
