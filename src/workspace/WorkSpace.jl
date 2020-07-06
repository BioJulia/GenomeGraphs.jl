
mutable struct WorkSpace
    sdg::Graphs.SequenceDistanceGraph{LongDNASeq}
    paired_reads_datastores::Vector{Union{PairedReads{DNAAlphabet{2}},PairedReads{DNAAlphabet{4}}}}
    long_reads_datastores::Vector{LongReads}
    linked_reads_datastores::Vector{LinkedReads}
    mer_count_stores::Vector{IndexedCounts}
end

"Create an empty workspace"
function WorkSpace()
    return WorkSpace(Graphs.SequenceDistanceGraph{LongDNASeq}(),
                     Vector{PairedReads}(),
                     Vector{LongReads}(),
                     Vector{LinkedReads}(),
                     Vector{IndexedCounts}())
end

graph(ws::WorkSpace) = ws.sdg

"""
    paired_reads(ws::WorkSpace, id::Symbol)

Get a reference to the `PairedReads` in the workspace that has the name
specified by `nm`.
"""
function paired_reads(ws::WorkSpace, nm::Symbol)
    for ds in ws.paired_reads_datastores
        if name(ds) == nm
            return ds
        end
    end
    error(string("WorkSpace has no paired read datastore called ", nm))
end

function add_paired_reads!(ws::WorkSpace, reads::PairedReads)
    push!(ws.paired_reads_datastores, reads)
    return ws
end

function add_paired_reads!(ws::WorkSpace, file::String, name::Union{String,Nothing} = nothing)
    reads = open(PairedReads, file, name)
    add_paired_reads!(ws, reads)
    return ws
end

function add_mer_counts!(ws::WorkSpace, ::Type{M}, mode::KmerCountMode, name::Symbol) where {M<:AbstractMer}
    push!(ws.mer_count_stores, MerCounts{M}(name, ws.sdg, mode))
    return ws
end

function mer_counts(ws::WorkSpace, name::Symbol)
    for mc in ws.mer_count_stores
        if mc.name == name
            return mc
        end
    end
    error(string("WorkSpace has no mer count datastore datastore called ", name))
end

function Base.show(io::IO, ws::WorkSpace)
    println("Graph Genome workspace")
    
    println(io, ' ', summary(ws.sdg))
    
    println(" Paired Read Datastores (", length(ws.paired_reads_datastores), "):")
    for ds in ws.paired_reads_datastores
        println(io, "  ", summary(ds)[23:end])
    end
    
    println(" Long Read Datastores (", length(ws.long_reads_datastores), "):")
    for ds in ws.long_reads_datastores
        println(io, "  ", summary(ds)[23:end])
    end
    
    println(" Linked Read Datastores (", length(ws.linked_reads_datastores), "):")
    for ds in ws.linked_reads_datastores
        println(io, "  ", summary(ds)[23:end])
    end
    
    println(" Mer Count Datastores (", length(ws.mer_count_stores), "):")
    for ds in ws.mer_count_stores
        msg = summary(ds)
        i = first(findfirst("'", msg))
        println(io, "  ", msg[i:end])
    end
end