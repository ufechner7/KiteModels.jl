using KiteModels
using Pkg
if ! ("Documenter" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end

using Documenter

DocMeta.setdocmeta!(KiteModels, :DocTestSetup, :(using KitePodSimulator); recursive=true)

makedocs(;
    modules=[KiteModels],
    authors="Uwe Fechner <fechner@aenarete.eu> and contributors",
    repo="https://github.com/ufechner7/KiteModels.jl/blob/{commit}{path}#{line}",
    sitename="KiteModels.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ufechner7.github.io/KiteModels.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Types" => "types.md",
        "Functions" => "functions.md",
        "Parameters" => "parameters.md",
        "Examples_1p" => "examples.md",
        "Examples_4p" => "examples_4p.md",
    ],
)

deploydocs(;
    repo="github.com/ufechner7/KiteModels.jl",
    devbranch="main",
)
