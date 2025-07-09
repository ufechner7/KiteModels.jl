using KiteModels, KiteUtils, Test

set = load_settings("system.yaml")
kps3::KPS3 = KPS3(set)
kps4::KPS4 = KPS4(set)
set = load_settings("system_ram.yaml")
saw::SymbolicAWEModel = SymbolicAWEModel(set)

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
    @test saw isa KiteUtils.AbstractKiteModel
    @test saw.set isa KiteUtils.Settings
    @test isnothing(saw.integrator)
    @test saw.am isa AtmosphericModel
    @test saw.iter isa Int64
end
nothing
