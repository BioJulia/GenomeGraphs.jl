using Documenter, BioSequenceGraphs

makedocs(
    format = Documenter.HTML(),
    sitename = "GenomeGraphs.jl",
    authors = "Ben J. Ward & Arda Akdemir",
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Guide" => "man/guide.md"
        ]
    ],
    
)

deploydocs(
    repo = "github.com/BioJulia/BioSequenceGraphs.jl.git",
    deps = nothing,
    make = nothing
)