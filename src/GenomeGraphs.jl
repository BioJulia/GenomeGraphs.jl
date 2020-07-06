__precompile__()

module GenomeGraphs

export
    ###
    ### Re-exports of GenomeGraphs framework sub-components.
    ###
    
    # Re-exports of ReadDatastores
    ReadDatastore,
    PairedReads,
    LongReads,
    LinkedReads,
    buffer,
    FwRv,
    
    # Re-exports of FASTX
    FASTA,
    FASTQ,
    
    ###
    ### Re-exports of BioSequences.jl
    ###
    DNAAlphabet,
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
    ### Re-exports of KmerAnalysis.jl
    ###
    mer,
    freq,
    Canonical,
    NonCanonical,
    CANONICAL,
    NONCANONICAL,
    # Kmer counters
    serial_mem,
    dist_mem,
    spectra,
    
    empty_graph,
    #GRAPH_TYPE,
    
    ### WorkSpace
    #WorkSpace,
    #add_paired_reads!,
    #paired_reads,
    #add_mer_counts!,
    #mer_counts,
    read_datastore,
    ###
    ### Processes
    ###
    dbg,
    dbg!,
    remove_tips!

using BioSequences, FASTX, ReadDatastores, KmerAnalysis
import BioSequences.EveryMerIterator

include("Graphs.jl")       # Submodule defining the key Graph type and basic methods.
include("GraphIndexes.jl") # Submodule defining types that allow indexing of a graph.



#include("workspace/WorkSpace.jl")
#include("views/NodeView.jl")

# Utility function for more quickly making a read datastore.
function read_datastore(R1file::String, R2file::String, name::String, minlen::Int, maxlen::Int, insertlen::Int, mode::PairedReadOrientation)
    fwq = open(FASTQ.Reader, R1file)
    rvq = open(FASTQ.Reader, R2file)
    return PairedReads{DNAAlphabet{4}}(fwq, rvq, name, name, minlen, maxlen, insertlen, mode)
end

empty_graph(::Type{T}) where {T<:LongSequence} = Graphs.SequenceDistanceGraph{T}()

include("dbg.jl")
include("remove_tips.jl")

end # module GenomeGraphs
