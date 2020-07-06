using Documenter, GenomeGraphs, Pkg

makedocs(
    modules = [GenomeGraphs, GenomeGraphs.Graphs, GenomeGraphs.GraphIndexes],
    format = Documenter.HTML(),
    sitename = "GenomeGraphs.jl",
    authors = replace(join(Pkg.TOML.parsefile("Project.toml")["authors"], ", "), r" <.*?>" => "" ),
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Installation" => "man/install.md"
            "Genome Assembly" => "man/assembly.md"
        ],
        "API" => [
            "Graphs submodule" => "api/Graphs.md",
            "GraphIndexes submodule" => "api/GraphIndexes.md"
        ]
    ],
    
)

deploydocs(
    repo = "github.com/BioJulia/GenomeGraphs.jl.git",
    push_preview = true,
    deps = nothing,
    make = nothing
)