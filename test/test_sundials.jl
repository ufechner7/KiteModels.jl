# SPDX-FileCopyrightText: 2025 Uwe Fechner
# SPDX-License-Identifier: MIT

# unittest for Sundials.jl, works fine with 4.11.4 and older
# also requires v6.115.4 of DiffEqBase or older
using StaticArrays, Sundials, Test

const KVec3    = MVector{3, Float64}

Base.@kwdef mutable struct KPS4{S, T, P}
    v_apparent::T =       zeros(S, 3)
end

const kps4 = KPS4{Float64, KVec3, 6+4+1}()

function residual!(res, yd, y::MVector{S, Float64}, s::KPS4, time) where S
end

function init!(t_end=1.0)
    y0 = MVector{62, Float64}([13.970413450119487, 0.0, 21.238692070636343, 27.65581376097752, 0.0, 42.66213714321849, 40.976226230518435, 0.0, 64.314401166278, 53.87184032029182, 0.0, 86.22231803750196, 66.28915240374937, 0.0, 108.4048292516046, 78.17713830204762, 0.0, 130.87545423106485, 79.56930502428155, 0.0, 135.70836376062155, 80.90383289255747, 0.0, 137.7696816741141, 80.60126812407692, 2.4016533873456325, 135.3023287520457, 80.60126812407692, -2.4016533873456325, 135.3023287520457, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 150.0, 0.0])
    yd0= MVector{62, Float64}([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -9.81, 0.0, 0.0, -9.81, 0.0, 0.0, -9.81, 0.0, 0.0, -9.81, 0.0, 0.0, -9.81, 0.0, 0.0, -9.81, 0.0, 0.0, -9.81, 0.0, 0.0, -9.81, 0.0, 0.0, -9.81, 0.0, 0.0, -9.81, 0.0, 0.0])
    differential_vars = ones(Bool, length(y0))
    solver  = IDA(linear_solver=Symbol("GMRES"), max_order = 4)
    tspan   = (0.0, t_end) 
    abstol  = 0.0006 # max error in m/s and m
    prob    = DAEProblem(residual!, yd0, y0, tspan, kps4, differential_vars=differential_vars)
    integrator = Sundials.init(prob, solver, abstol=abstol, reltol=0.001)
end

@testset "test_simulate        " begin
    STEPS = 500
    integrator = init!()
end
nothing
