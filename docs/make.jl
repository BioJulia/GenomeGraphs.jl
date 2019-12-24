using Documenter, GenomeGraphs, Pkg

makedocs(
    modules = [GenomeGraphs, GenomeGraphs.Graphs, GenomeGraphs.MerTools],
    format = Documenter.HTML(),
    sitename = "GenomeGraphs.jl",
    authors = replace(join(Pkg.TOML.parsefile("Project.toml")["authors"], ", "), r" <.*?>" => "" ),
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Guide" => "man/guide.md"
        ],
        "API" => [
            "Graphs submodule" => "api/Graphs.md",
            "MerTools submodule" => "api/MerTools.md"
        ]
    ],
    
)

deploydocs(
    repo = "github.com/BioJulia/GenomeGraphs.jl.git",
    push_preview = true,
    deps = nothing,
    make = nothing
)