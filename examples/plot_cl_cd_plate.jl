# plot the lift and drag coefficients as function of angle of attack
# of any of the plates of the kite

using KiteModels

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
end
plot(ALPHA, CL)