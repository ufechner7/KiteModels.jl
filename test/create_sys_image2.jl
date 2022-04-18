using Pkg
if ! ("PackageCompiler" ∈ keys(Pkg.project().dependencies))
    @info "Installing PackageCompiler ..."
    Pkg.add("PackageCompiler")
end
if ! ("Plots" ∈ keys(Pkg.project().dependencies))
    @info "Installing Plots ..."
    Pkg.add("Plots")
end
@info "Loading packages ..."
using KiteUtils, KitePodModels, KiteModels, Plots
using PackageCompiler

@info "Creating sysimage ..."
push!(LOAD_PATH,joinpath(pwd(),"src"))

PackageCompiler.create_sysimage(
    [:KiteUtils, :KitePodModels, :KiteModels, :Plots];
    sysimage_path="kps-image_tmp.so",
    precompile_execution_file=joinpath("test", "test_for_precompile.jl")
)