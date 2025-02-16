# plot the lift and drag coefficients as function of angle of attack
# of any of the plates of the kite

using KiteModels

if haskey(ENV, "USE_V9")
    set = deepcopy(load_settings("system_v9.yaml"))
else
    set = deepcopy(load_settings("system.yaml"))
end

using Pkg
if ! ("Test" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using ControlPlots, LaTeXStrings
plt.close("all")

set.v_wind = 14 # 25
kcu::KCU = KCU(set)
kps4::KPS4 = KPS4(kcu)

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

ALPHA = -10:0.1:20
plot_cl_cd(ALPHA)

cl1 = kps4.calc_cl(-5)
cl2 = kps4.calc_cl(12)
dcl_over_dalpha= (cl2-cl1)/deg2rad(12+5)
println("dCL/dalpha = ", round(dcl_over_dalpha, digits=3))