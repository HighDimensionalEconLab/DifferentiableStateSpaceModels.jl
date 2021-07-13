using DifferentiableStateSpaceModels
using Documenter

DocMeta.setdocmeta!(
    DifferentiableStateSpaceModels,
    :DocTestSetup,
    :(using DifferentiableStateSpaceModels);
    recursive = true,
)

makedocs(;
    modules = [DifferentiableStateSpaceModels],
    authors = "Jesse Perla <HighDimensionalEconLab.com> and contributors",
    repo = "https://github.com/HighDimensionalEconLab/DifferentiableStateSpaceModels.jl/blob/{commit}{path}#{line}",
    sitename = "DifferentiableStateSpaceModels.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://HighDimensionalEconLab.github.io/DifferentiableStateSpaceModels.jl",
        assets = String[],
    ),
    pages = ["Home" => "index.md"],
)

deploydocs(; repo = "github.com/HighDimensionalEconLab/DifferentiableStateSpaceModels.jl")
