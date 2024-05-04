using Pkg
if ! ("Plots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using Plots
Plots.__init__()
using KiteModels, KitePodModels, KiteUtils

set = deepcopy(se())
kcu::KCU = KCU(set)
kps4::KPS4 = KPS4(kcu)

include("../examples/plot2d.jl")

reltime = 0.0
integrator = KiteModels.init_sim!(kps4, stiffness_factor=0.5)
plot2d(kps4.pos, reltime; zoom=false, front=false, segments=set.segments)