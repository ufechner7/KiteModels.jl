using Pkg
if ! ("PackageCompiler" ∈ keys(Pkg.project().dependencies))
    @info "Installing PackageCompiler ..."
    Pkg.add("PackageCompiler")
end
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    @info "Installing ControlPlots ..."
    Pkg.add("ControlPlots")
end
@info "Loading packages ..."
using KiteUtils, KitePodModels, KiteModels, ControlPlots, AtmosphericModels, Reexport
using PackageCompiler

@info "Creating sysimage ..."
push!(LOAD_PATH,joinpath(pwd(),"src"))

PackageCompiler.create_sysimage(
    [:KiteUtils, :KitePodModels, :KiteModels, :ControlPlots, :AtmosphericModels, :Reexport];
    sysimage_path="kps-image_tmp.so",
    precompile_execution_file=joinpath("test", "test_for_precompile.jl")
)