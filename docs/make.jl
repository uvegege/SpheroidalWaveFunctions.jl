using Documenter
using SpheroidalWaveFunctions

makedocs(modules = [SpheroidalWaveFunctions],
    sitename = "SpheroidalWaveFunctions.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md"
    ]
)

deploydocs(
    repo = "github.com/uvegege/SpheroidalWaveFunctions.jl.git",
    devbranch = "main"
)
