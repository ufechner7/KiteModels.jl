# SPDX-FileCopyrightText: 2025 Uwe Fechner
#
# SPDX-License-Identifier: MIT

using KiteModels
using Test

set_data_path(joinpath(dirname(dirname(pathof(KiteModels))), "data"))
set = deepcopy(load_settings("system.yaml"))

TEST_KCU = KCU(set)
S = KPS4(TEST_KCU)

@testset verbose=true "test_rotational_inertia" begin
    @testset "kite including KCU around CoM " begin
        Ixx, Ixy, Ixz, Iyy, Iyz, Izz = calculate_rotational_inertia!(S, true, false)

        @test isapprox(Ixx, 124.53,     rtol=0.005)
        @test isapprox(Ixy, 0,          rtol=0.001)
        @test isapprox(Ixz, -8.805,     rtol=0.007)
        @test isapprox(Iyy, 111.227,    rtol=0.006)
        @test isapprox(Iyz, 0,          rtol=0.001)
        @test isapprox(Izz, 19.516,     rtol=0.001)
    end

    @testset "kite no KCU around CoM" begin
        Ixx, Ixy, Ixz, Iyy, Iyz, Izz = calculate_rotational_inertia!(S, false, false)

        @test isapprox(Ixx, 21.5607,    rtol=0.001)
        @test isapprox(Ixy, 0,          rtol=0.001)
        @test isapprox(Ixz, 1.589,      rtol=0.001)
        @test isapprox(Iyy, 7.2073,     rtol=0.001)
        @test isapprox(Iyz, 0,          rtol=0.001)
        @test isapprox(Izz, 18.4558,    rtol=0.001)
    end

    @testset "kite including KCU around KCU" begin
        Ixx, Ixy, Ixz, Iyy, Iyz, Izz = calculate_rotational_inertia!(S, true, true)

        @test isapprox(Ixx, 200.533,    rtol=0.001) 
        @test isapprox(Ixy, 0,          rtol=0.001) 
        @test isapprox(Ixz, -16.478,    rtol=0.001)
        @test isapprox(Iyy, 188.004,    rtol=0.001)
        @test isapprox(Iyz, 0,          rtol=0.001)
        @test isapprox(Izz, 20.2907,    rtol=0.001)
    end

    @testset "kite no KCU around KCU" begin
        Ixx, Ixy, Ixz, Iyy, Iyz, Izz = calculate_rotational_inertia!(S, false, true)

        @test isapprox(Ixx, 200.533,    rtol=0.001) 
        @test isapprox(Ixy, 0,          rtol=0.001) 
        @test isapprox(Ixz, -16.478,    rtol=0.001)
        @test isapprox(Iyy, 188.004,    rtol=0.001)
        @test isapprox(Iyz, 0,          rtol=0.001)
        @test isapprox(Izz, 20.2907,    rtol=0.001)
    end
end