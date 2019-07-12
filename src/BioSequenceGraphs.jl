__precompile__()

module BioSequenceGraphs

export
    ### Sequence Distance Graph
    SequenceDistanceGraph,
    # Basic queries and properties
    nodes,
    n_nodes,
    each_node_id,
    node,
    links,
    sequence,
    # Graph traversal
    get_next_nodes,
    get_previous_nodes,
    get_all_unitigs

using BioSequences

include("graph/SequenceDistanceGraph.jl")
include("graph/graph_building.jl")
end # module BioSequenceGraphs
