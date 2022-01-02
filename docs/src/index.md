```@meta
CurrentModule = KiteModels
```

# KiteModels

Documentation for the package [KiteModels](https://github.com/ufechner7/KiteModels.jl).

## Background


## Installation

Download [Julia 1.6](http://www.julialang.org) or later, if you haven't already. You can add KiteModels from  Julia's package manager, by typing 
```
] add KiteModels
``` 
at the Julia prompt.

If you are using Windows, it is suggested to install git and bash, too. This is explained for example here: [Julia on Windows](https://github.com/ufechner7/KiteViewer/blob/main/doc/Windows.md) .

## Provides

The type KPS3 with the **residual!** function for a DAE solver, representing the model. 

Additional functions to provide inputs and outputs of the model on each time step. Per time step the residual! function is called as many times as needed to find the solution at the end
of the time step.


Click on **Functions** on the left to see the exported functions.

Author: Uwe Fechner (uwe.fechner.msc@gmail.com)
