

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
    plt.figure("spectrum")
    f_min = 1.5 # Hz
    f_max = 4.0 # Hz
    f_ex = []
    aoa_eff = []
    for i in 1:length(spectrum.f_ex)
        if spectrum.f_ex[i] < f_min
            continue
        end
        if spectrum.f_ex[i] > f_max
            break
        end
        push!(f_ex, spectrum.f_ex[i])
        push!(aoa_eff, spectrum.aoa_eff[i])
    end
    plt.plot(f_ex, todb.(aoa_eff); label=name)
    plt.xlabel("f_ex [Hz]")
    plt.ylabel("AOA amplitude [dB°]")
    plt.gca().set_xscale("log")
    plt.grid(true)
    plt.legend()
end

plot_spectrum("spectrum2_8.0_0.09")
plot_spectrum("spectrum2_8.0_0.0")