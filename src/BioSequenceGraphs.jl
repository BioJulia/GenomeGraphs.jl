__precompile__()

module BioSequenceGraphs

export
    ###
    ### Re-exports of GenomeGraphs framework sub-components.
    ###
    
    # ReadDatastores
    ReadDatastore,
    PairedReads,
    LongReads,
    LinkedReads,
    buffer,
    
    # FASTX
    FASTA,
    FASTQ,
    
    # BioSequences
    BioSequence,
    AbstractMer,
    Mer,
    BigMer,
    DNAMer,
    DNAKmer,
    

    ###
    ### Sequence Distance Graph
    ###
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
    
    ###
    ### MerFreq
    ###
    collapse_sorted!,
    collapse!,
    collapse_into_freqs!,
    collapse_into_freqs,
    collapse_into_freqs_sorted!,
    collapse_into_freqs_sorted,
    merge_into!,
    merge_into_sorted!,
    hist,
    hist!,
    
    
    
    ### WorkSpace
    WorkSpace,
    add_paired_reads!

using BioSequences, FASTX, ReadDatastores
import BioSequences.EveryMerIterator

include("mertools/MerFreq.jl")
include("mertools/counting.jl")

include("graph/SequenceDistanceGraph.jl")
include("graph/graph_building.jl")
include("workspace/WorkSpace.jl")
end # module BioSequenceGraphs
