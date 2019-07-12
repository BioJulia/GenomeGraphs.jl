using Documenter, BioSequenceGraphs

makedocs(
    format = Documenter.HTML(),
    sitename = "GenomeGraphs.jl",
    pages = [
        "Home" => "index.md",
        "Sequence Distance Graphs" => "graphs.md"
    ],
    authors = "Ben J. Ward & Arda Akdemir"
)