# activate the test environment if needed
using Pkg
if ! ("PackageCompiler" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
@info "Loading packages ..."
using Dierckx, StaticArrays, LinearAlgebra, Parameters, NLsolve, DocStringExtensions, Sundials, KiteUtils, KitePodModels
using PackageCompiler

@info "Creating sysimage ..."
push!(LOAD_PATH,joinpath(pwd(),"src"))

PackageCompiler.create_sysimage(
    [:Dierckx, :StaticArrays, :Parameters, :NLsolve, :DocStringExtensions, :Sundials, :KiteUtils, :KitePodModels];
    sysimage_path="kps-image_tmp.so",
    precompile_execution_file=joinpath("test", "test_for_precompile.jl")
)