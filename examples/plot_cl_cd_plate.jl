# plot the lift and drag coefficients as function of angle of attack
# of any of the plates of the kite

using KiteModels

set = deepcopy(load_settings("system.yaml"))

using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using ControlPlots, LaTeXStrings

set.v_wind = 14 # 25
kcu = KCU(set)
kps4 = KPS4(kcu)

function plot_cl_cd(alpha)
    cl = zeros(length(alpha))
    cd = zeros(length(alpha))
    for (i, alpha) in pairs(ALPHA)
        cl[i] = kps4.calc_cl(alpha)
        cd[i] = kps4.calc_cd(alpha)
    end
    display(plot(ALPHA, [cl, cd]; xlabel=L"\mathrm{AoA}~\alpha", ylabel="CL, CD", labels=["CL", "CD"], fig="CL_CD"))
    display(plot(ALPHA, [cl./cd]; xlabel=L"\mathrm{AoA}~\alpha", ylabel="LoD", fig="LoD"))
end

ALPHA = 0:0.1:14
plot_cl_cd(ALPHA)


