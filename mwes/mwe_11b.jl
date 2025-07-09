# SPDX-FileCopyrightText: 2025 Uwe Fechner
#
# SPDX-License-Identifier: MIT

using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using ControlPlots
using KiteModels, KitePodModels, KiteUtils

set = deepcopy(se())
kcu::KCU = KCU(set)
kps4::KPS4 = KPS4(kcu)

reltime = 0.0
integrator = KiteModels.init!(kps4, stiffness_factor=0.5)
plot2d(kps4.pos, reltime; zoom=true, front=false, segments=set.segments)