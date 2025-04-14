using KiteModels
using Timers
using Pkg 
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using KiteUtils
using Dierckx
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkit: Symbolics, @register_symbolic
using OrdinaryDiffEqCore
using ModelingToolkit, LinearAlgebra, Timers, Parameters, ControlPlots
using Pkg

set = load_settings("system_kps5.yaml")
