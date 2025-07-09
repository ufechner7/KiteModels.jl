using KiteModels, KiteUtils, Test

set = load_settings("system.yaml")
kps3::KPS3 = KPS3(set)
# kps4::KPS4 = KPS4(set)
# saw::SymbolicAWEModel = 

@testset "KPS3 constructor interface" begin
    @test kps3 isa KiteUtils.AbstractKiteModel
    @test kps3.set isa KiteUtils.Settings
    @test isnothing(kps3.integrator)
end
# @testset "KPS4 constructor interface" begin
#     @test kps4 isa KiteUtils.AbstractKiteModel
#     @test kps4.set isa KiteUtils.Settings
#     @test isnothing(kps4.integrator)
# end
nothing
