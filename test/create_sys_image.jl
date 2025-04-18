# activate the test environment if needed
using Pkg
if ! ("PackageCompiler" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
    Pkg.update()
end
@info "Loading packages ..."
using Dierckx, StaticArrays, LinearAlgebra, Parameters, NLsolve, DocStringExtensions, Sundials, KiteUtils, 
      KitePodModels, AtmosphericModels, OrdinaryDiffEqCore, OrdinaryDiffEqBDF, OrdinaryDiffEqSDIRK, ModelingToolkit,
      DSP, JLD2, Colors, REPL, VortexStepMethod, NonlinearSolve, FiniteDiff, OrdinaryDiffEqNonlinearSolve, StatsBase
using PackageCompiler, BenchmarkTools, Documenter

@info "Creating sysimage ..."
push!(LOAD_PATH,joinpath(pwd(),"src"))

PackageCompiler.create_sysimage(
    [:Dierckx, :StaticArrays, :Parameters, :NLsolve, :DocStringExtensions, :Sundials, :KiteUtils, 
     :KitePodModels, :AtmosphericModels, :OrdinaryDiffEqCore, :OrdinaryDiffEqBDF, :OrdinaryDiffEqSDIRK, 
     :OrdinaryDiffEqNonlinearSolve, :StatsBase, :PackageCompiler, :BenchmarkTools, :Documenter,
     :ModelingToolkit, :JLD2, :Colors, :REPL, :VortexStepMethod, :NonlinearSolve, :DSP, :FiniteDiff];
    sysimage_path="kps-image_tmp.so",
    include_transitive_dependencies=true,
    precompile_execution_file=joinpath("test", "test_for_precompile.jl")
)