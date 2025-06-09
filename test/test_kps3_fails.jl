# SPDX-FileCopyrightText: 2025 Uwe Fechner
# SPDX-License-Identifier: MIT

using Pkg
if ! ("Test" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end

using Test, BenchmarkTools, StaticArrays, LinearAlgebra, KiteUtils

using KiteModels, KitePodModels

const SEGMENTS = se().segments
if ! @isdefined kcu
    const kcu = KCU(se())
    const kps = KPS3(kcu)
end
res1 = zeros(SVector{SEGMENTS+1, KiteModels.KVec3})
res2 = deepcopy(res1)
if ! @isdefined res3
    const res3 = vcat(reduce(vcat, vcat(res1, res2)), zeros(2))
end

@testset verbose = true "KPS3 tests...." begin

function set_defaults()
    KiteModels.clear!(kps)
    kps.set.l_tethers[1] = 150.0
    kps.set.elevation = 60.0
    kps.set.area = 20.0
    kps.set.rel_side_area = 50.0
    kps.set.v_wind = 8.0
    kps.set.mass = 11.4
    kps.set.damping =  2 * 473.0
    kps.set.alpha = 1.0/7
    kps.set.c_s = 0.6
end

set_defaults()

@testset "calc_cl              " begin
    @test isapprox(kps.calc_cl(-5.0), 0.150002588978, atol=1e-4) 
    @test isapprox(kps.calc_cl( 0.0), 0.200085035326, atol=1e-4) # fails
    @test isapprox(kps.calc_cl(10.0), 0.574103590856, atol=1e-4) # fails
    @test isapprox(kps.calc_cl(20.0), 1.0, atol=1e-4)
end
end

println("alpha_cl: $(se().alpha_cl), cl_list: $(se().cl_list)")
nothing
