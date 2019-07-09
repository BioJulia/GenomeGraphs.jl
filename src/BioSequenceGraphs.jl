__precompile__()

module BioSequenceGraphs

export SequenceDistanceGraph,
    new_graph_from_kmerlist

using BioSequences

include("graph/SequenceDistanceGraph.jl")
include("graph/graph_building.jl")
end # module BioSequenceGraphs
