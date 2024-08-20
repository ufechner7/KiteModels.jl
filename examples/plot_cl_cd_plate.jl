# plot the lift and drag coefficients as function of angle of attack
# of any of the plates of the kite

using KiteModels, LaTeXStrings

set = deepcopy(se())

using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using ControlPlots

plot(set.alpha_cl, set.cl_list)

ALPHA = -180:1:180
CL = zeros(length(ALPHA))
CD = zeros(length(ALPHA))
for (i, alpha) in pairs(ALPHA)
    CL[i] = KiteModels.calc_cl(alpha)
    CD[i] = KiteModels.calc_cd(alpha)
end
display(plot(ALPHA, [CL, CD]; xlabel=L"\mathrm{AoA}~\alpha", ylabel="CL, CD", labels=["CL", "CD"]))