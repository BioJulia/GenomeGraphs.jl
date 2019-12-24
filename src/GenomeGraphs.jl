__precompile__()

module GenomeGraphs

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
    FwRv,
    
    # FASTX
    FASTA,
    FASTQ,
    
    # BioSequences
    BioSequence,
    LongSequence,
    AbstractMer,
    Mer,
    BigMer,
    DNAMer,
    DNAKmer,
    BigDNAMer,
    BigDNAKmer,
    
    ###
    ### MerCounts
    ###
    MerCounts,
    
    ### WorkSpace
    WorkSpace,
    add_paired_reads!,
    paired_reads,
    add_mer_counts!,
    mer_counts,
    
    ###
    ### Processes
    ###
    dbg,
    dbg!,
    remove_tips!

include("MerTools.jl")   # MerTools submodule.
include("Graphs.jl")     # Graphs submodule.

using BioSequences, FASTX, ReadDatastores
import BioSequences.EveryMerIterator

include("indexes/unique-kmers.jl")

include("datastores/kmer-counts.jl")
include("workspace/WorkSpace.jl")
include("views/NodeView.jl")

include("processes/dbg.jl")
include("processes/remove_tips.jl")
end # module GenomeGraphs
