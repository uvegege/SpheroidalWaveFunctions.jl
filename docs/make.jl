using Documenter
using SpheroidalWaveFunctions

makedocs(
    sitename = "SpheroidalWaveFunctions.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    modules = [SpheroidalWaveFunctions],
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md"
    ]
)

deploydocs(
    repo = "github.com/uvegege/SpheroidalWaveFunctions.jl.git",
    devbranch = "main"
)
