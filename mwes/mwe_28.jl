using KiteModels, KiteUtils, Test

set = load_settings("system.yaml")
kcu::KCU  = KCU(set)
kps3::KPS3 = KPS3(kcu)

@testset "KPS3 constructor interface" begin
    @test kps3 isa KiteUtils.AbstractKiteModel
    @test kps3.set isa KiteUtils.Settings
    # kps3.integrator
end
nothing
