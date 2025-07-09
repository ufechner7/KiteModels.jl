# Copyright (c) 2022, 2024 Uwe Fechner
# SPDX-License-Identifier: MIT

using KiteModels, KitePodModels

set_data_path("data")
set = load_settings("system_v9.yaml")
set.elevation = 70.8

kcu::KCU = KCU(set)
kps4::KPS4 = KPS4(kcu)

clear!(kps4)
KiteModels.set_depower_steering!(kps4, kps4.set.depower_offset/100.0, 0.0)
@time KiteModels.init!(kps4; delta=0.001, stiffness_factor=0.5, prn=false)

println("\nlift, drag    [N]  : $(KiteModels.lift_drag(kps4))")
println("\nSpring forces:")
spring_forces(kps4)