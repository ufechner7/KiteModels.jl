using Pkg
if ! ("Plots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using Plots
Plots.__init__()
using KiteModels, KitePodModels, KiteUtils

set = deepcopy(se())
kcu::KCU = KCU(set)
kps3::KPS3 = KPS3(kcu)

include("../examples/plot2d.jl")

reltime = 0.0
integrator = KiteModels.init_sim!(kps3, stiffness_factor=0.04)
plot2d(kps3.pos, reltime; zoom=false, front=false, segments=set.segments)
# display(p)   