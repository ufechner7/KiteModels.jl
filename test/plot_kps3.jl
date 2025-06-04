# Copyright (c) 2022, 2024, 2025 Uwe Fechner
# SPDX-FileCopyrightText: 2025 Uwe Fechner
#
# SPDX-License-Identifier: MIT

# activate the test environment if needed
using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end

using KiteModels
using KitePodModels
using KiteUtils
using ControlPlots
using LinearAlgebra

kcu::KCU = KCU(se())
kps3::KPS3 = KPS3(kcu)

FRONT_VIEW = false

clear!(kps3)
kps3.stiffness_factor = 0.04

y0, yd0 = KiteModels.init(kps3)

find_steady_state!(kps3, prn=false)

println("kite distance: $(norm(kps3.pos[end]))")
println(KiteModels.spring_forces(kps3))
println("alpha_depower [deg]: $(rad2deg(kps3.alpha_depower))")
println("lift, drag    [N]  : $(KiteModels.lift_drag(kps3))")

plot2d(kps3.pos; zoom=false, front=FRONT_VIEW)