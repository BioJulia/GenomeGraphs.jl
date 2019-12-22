using Documenter, GenomeGraphs

makedocs(
    modules = [GenomeGraphs, GenomeGraphs.Graphs],
    format = Documenter.HTML(),
    sitename = "GenomeGraphs.jl",
    authors = "Ben J. Ward & Arda Akdemir",
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Guide" => "man/guide.md"
        ],
        "API" => [
            "Graphs submodule" => "api/Graphs.md"
        ]
    ],
    
)

deploydocs(
    repo = "github.com/BioJulia/GenomeGraphs.jl.git",
    deps = nothing,
    make = nothing
)