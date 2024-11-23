

using Printf
using KiteModels, StatsBase, LinearAlgebra, DSP

if haskey(ENV, "USE_V9")
    set = deepcopy(load_settings("system_v9.yaml"))
else
    set = deepcopy(load_settings("system.yaml"))
end

using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using ControlPlots, JLD2
plt.close("all")

function plot_spectrum(name)
    todb(mag) = 20 * log10(mag)
    spectrum = jldopen("data/" * name * ".jld2") do file
        read(file, "spectrum")
    end
    plt.figure(name)
    plt.plot(spectrum.f_ex, todb.(spectrum.aoa_eff))
    plt.xlabel("f_ex [Hz]")
    plt.ylabel("AOA amplitude [dB°]")
    plt.gca().set_xscale("log")
    plt.grid(true)
end

plot_spectrum("spectrum_8.0_0.2")