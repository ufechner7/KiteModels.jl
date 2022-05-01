# KiteModels

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ufechner7.github.io/KiteModels.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ufechner7.github.io/KiteModels.jl/dev)
[![Build Status](https://github.com/ufechner7/KiteModels.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ufechner7/KiteModels.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/ufechner7/KiteModels.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ufechner7/KiteModels.jl)

## Kite power system models, consisting of tether and kite

## One point model
This model assumes the kite to be a point mass. This is sufficient to model the aerodynamic forces, but the dynamic with respect to turning action of the kite is not realistic.
When combined with an controller for the turn rate it can be used to simulate a pumping kite power system with medium accuracy.

## Four point model
This model assumes the kite to consist of four point masses. It reacts much more realistically than the one point model.

## Tether
The tether is modeled as point masses, connected by spring-damper elements. Aerodynamic drag is modelled realistically. When reeling out or in the unstreched length of the spring-damper elements
is varied. This does not translate into physics directly, but it avoids adding point masses at run-time, which would be even worse because it would introduce discontinuities. When using
Dyneema or simular high strenght materials for the tether the resulting system is very stiff which is a challenge for the solver.

## Further reading
These models are described in detail in [Dynamic Model of a Pumping Kite Power System](http://arxiv.org/abs/1406.6218).

## See also
- [Research Fechner](https://research.tudelft.nl/en/publications/?search=Uwe+Fechner&pageSize=50&ordering=rating&descending=true) for the scientic background of this code
- The application [KiteViewer](https://github.com/ufechner7/KiteViewer)
- the package [KiteUtils](https://github.com/ufechner7/KiteUtils.jl)
- the package [KitePodModels](https://github.com/aenarete/KitePodModels.jl)
- the package [KiteControllers](https://github.com/aenarete/KiteControllers.jl)

**Documentation** [Stable Version](https://ufechner7.github.io/KiteModels.jl/stable) --- [Development Version](https://ufechner7.github.io/KiteModels.jl/dev)


Author: Uwe Fechner (uwe.fechner.msc@gmail.com)