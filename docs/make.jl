# SPDX-FileCopyrightText: 2025 Uwe Fechner, Bart van de Lint
# SPDX-License-Identifier: MIT

using Pkg
if ("TestEnv" ∈ keys(Pkg.project().dependencies))
    if ! ("Documents" ∈ keys(Pkg.project().dependencies))
        using TestEnv; TestEnv.activate()
    end
end
using ControlPlots, VortexStepMethod
using KiteModels
using Documenter

DocMeta.setdocmeta!(KiteModels, :DocTestSetup, :(using KiteModels); recursive=true)

makedocs(;
    modules=[KiteModels],
    authors="Uwe Fechner <fechner@aenarete.eu>, Bart van de Lint <bart@vandelint.net> and contributors",
    repo="https://github.com/OpenSourceAWE/KiteModels.jl/blob/{commit}{path}#{line}",
    sitename="KiteModels.jl",
    format=Documenter.HTML(;
        repolink = "https://github.com/OpenSourceAWE/KiteModels.jl",
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://OpenSourceAWE.github.io/KiteModels.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Types" => "types.md",
        "Functions" => "functions.md",
        "SymbolicAWEModel" => "ram_air_kite.md",
        "Parameters" => "parameters.md",
        "Examples 1p" => "examples.md",
        "Examples 4p" => "examples_4p.md",
        "Examples SymbolicAWEModel" => "examples_ram_air.md",
        "SystemStructure for custom models" => "tutorial_system_structure.md",
        "Quickstart" => "quickstart.md",
        "Advanced usage" => "advanced.md",
    ],
)

deploydocs(;
    repo="github.com/OpenSourceAWE/KiteModels.jl",
    devbranch="main",
)
