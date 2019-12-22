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

include("Graphs.jl")

using BioSequences, FASTX, ReadDatastores
import BioSequences.EveryMerIterator

include("mertools/MerFreq.jl")
include("mertools/counting.jl")

include("indexes/unique-kmers.jl")

include("datastores/kmer-counts.jl")
include("workspace/WorkSpace.jl")
include("views/NodeView.jl")

include("processes/dbg.jl")
include("processes/remove_tips.jl")
end # module GenomeGraphs
