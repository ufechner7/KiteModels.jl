# SPDX-FileCopyrightText: 2025 Uwe Fechner
# SPDX-License-Identifier: MIT

using KiteModels, KiteUtils, Test

# set_data_path("data") 
set = deepcopy(load_settings("system.yaml"))
kps3::KPS3 = KPS3(set)
kps4::KPS4 = KPS4(set)
set_ram = deepcopy(load_settings("system_ram.yaml"))
sam::SymbolicAWEModel = SymbolicAWEModel(set_ram)

@testset "KPS3 constructor interface" begin
    @test kps3 isa KiteUtils.AbstractKiteModel
    @test kps3.set isa KiteUtils.Settings
    @test isnothing(kps3.integrator)
    @test kps3.am isa AtmosphericModel
    @test kps3.iter isa Int64
end
@testset "KPS4 constructor interface" begin
    @test kps4 isa KiteUtils.AbstractKiteModel
    @test kps4.set isa KiteUtils.Settings
    @test isnothing(kps4.integrator)
    @test kps4.am isa AtmosphericModel
    @test kps4.iter isa Int64
end
@testset "SymbolicAWEModel constructor interface" begin
    @test sam isa KiteUtils.AbstractKiteModel
    @test sam.set isa KiteUtils.Settings
    @test isnothing(sam.integrator)
    @test sam.am isa AtmosphericModel
    @test sam.iter isa Int64
end

@testset "KPS3 init_sim! interface" begin
    kps4.set.upwind_dir = rad2deg(-pi/2)
    integ = init_sim!(kps3; stiffness_factor=0.5, delta=0.0001, prn=false)
    @test integ isa KiteModels.OrdinaryDiffEqCore.ODEIntegrator
    @test kps3.integrator isa KiteModels.OrdinaryDiffEqCore.ODEIntegrator
end

@testset "KPS4 init_sim! interface" begin
    kps4.set.upwind_dir = rad2deg(-pi/2)
    integ = init_sim!(kps4; stiffness_factor=0.5, delta=0.0001, prn=false)
    @test integ isa KiteModels.OrdinaryDiffEqCore.ODEIntegrator
    @test kps4.integrator isa KiteModels.OrdinaryDiffEqCore.ODEIntegrator
end

@testset "SymbolicAWEModel init_sim! interface" begin
    sam.set.upwind_dir = rad2deg(-pi/2)
    integ = init_sim!(sam; stiffness_factor=0.5, delta=0.0001, prn=false)
    @test integ isa KiteModels.OrdinaryDiffEqCore.ODEIntegrator
    @test kps4.integrator isa KiteModels.OrdinaryDiffEqCore.ODEIntegrator
end
nothing