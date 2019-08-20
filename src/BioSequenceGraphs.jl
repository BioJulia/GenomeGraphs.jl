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
    get_all_unitigs,
    ### WorkSpace
    WorkSpace,
    add_paired_reads!

using BioSequences, FASTX, ReadDatastores

include("graph/SequenceDistanceGraph.jl")
include("graph/graph_building.jl")
include("workspace/WorkSpace.jl")
end # module BioSequenceGraphs
