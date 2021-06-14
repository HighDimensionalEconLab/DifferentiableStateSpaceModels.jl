using DifferentiableStateSpaceModels
using Documenter

DocMeta.setdocmeta!(DifferentiableStateSpaceModels, :DocTestSetup, :(using DifferentiableStateSpaceModels); recursive=true)

makedocs(;
    modules=[DifferentiableStateSpaceModels],
    authors="Jesse Perla <jesseperla@gmail.com> and contributors",
    repo="https://github.com/jlperla/DifferentiableStateSpaceModels.jl/blob/{commit}{path}#{line}",
    sitename="DifferentiableStateSpaceModels.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://jlperla.github.io/DifferentiableStateSpaceModels.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jlperla/DifferentiableStateSpaceModels.jl",
)
