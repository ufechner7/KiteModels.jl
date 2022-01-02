# KiteModels

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ufechner7.github.io/KiteModels.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ufechner7.github.io/KiteModels.jl/dev)
[![Build Status](https://github.com/ufechner7/KitePodSimulator.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ufechner7/KiteModels.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/ufechner7/KitePodSimulator.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ufechner7/KiteModels.jl)

## Kite power system models, consisting of tether and kite

## One point model
This model assumes the kite to be a point mass. This is sufficient to model the aerodynamic forces, but the dynamic with respect to turning action of the kite is not realistic.
When combined with an controller for the turn rate it can be used to simulate a pumping kite power system with medium accuracy.

## Four point model
This model assumes the kite to consist of four point masses. It reacts much more realistically than the one point model.

## Tether
