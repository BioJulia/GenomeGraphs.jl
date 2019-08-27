
mutable struct WorkSpace
    sdg::SequenceDistanceGraph{LongDNASeq}
    paired_reads_datastores::Vector{PairedReads}
    long_reads_datastores::Vector{LongReads}
    linked_reads_datastores::Vector{LinkedReads}
    #mer_count_stores::Vector{MerCounts}
end

"Create an empty workspace"
function WorkSpace()
    return WorkSpace(SequenceDistanceGraph{LongDNASeq}(),
                     Vector{PairedReads}(),
                     Vector{LongReads}(),
                     Vector{LinkedReads}())
end

"""
    paired_reads(ws::WorkSpace, id::String)

Get a reference to the `PairedReads` in the workspace that has the name
specified by `nm`.
"""
function paired_reads(ws::WorkSpace, nm::String)
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
end