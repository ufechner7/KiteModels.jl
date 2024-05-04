using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
    # pkg"add ControlPlots#main"
end
using ControlPlots
using KiteModels, KitePodModels, KiteUtils

set = deepcopy(se())
kcu::KCU = KCU(set)
kps4::KPS4 = KPS4(kcu)

include("../examples/plot2D.jl")

reltime = 0.0
integrator = KiteModels.init_sim!(kps4, stiffness_factor=0.5)

lines, sc, txt = nothing, nothing, nothing
plot2d(kps4.pos, reltime; zoom=true, front=false, segments=set.segments, lines, sc, txt)
nothing