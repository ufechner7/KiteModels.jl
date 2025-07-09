using KiteModels, Test

set = load_settings("system.yaml")
kcu::KCU  = KCU(set)
kps3::KPS3 = KPS3(kcu)

@test kps3 isa AbstractKiteModel
