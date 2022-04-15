# activate the test environment if needed
using Pkg
if ! ("Plots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end

using KiteModels
using KitePodModels
using KiteUtils
using Plots
using LinearAlgebra

const kcu = KCU()
const kps4 = KPS3(kcu)

include("plot_initial_state.jl")