# SPDX-FileCopyrightText: 2025 Uwe Fechner
#
# SPDX-License-Identifier: MIT

using StaticArrays, LinearAlgebra, KiteUtils
using KiteModels, KitePodModels

set_data_path(joinpath(dirname(dirname(pathof(KiteModels))), "data"))
set = deepcopy(load_settings("system.yaml"))
kcu::KCU = KCU(set)

kps4::KPS4 = KPS4(kcu)

dt = 0.05

clear!(kps)
KiteModels.set_depower_steering!(kps, kps.set.depower_offset/100.0, 0.0)
kps.stiffness_factor = 0.5

@time KiteModels.find_steady_state!(kps, prn=true)

println("\nlift, drag    [N]  : $(KiteModels.lift_drag(kps))")
# println("\nSpring forces:")
# spring_forces(kps)
