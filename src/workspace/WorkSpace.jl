
mutable struct WorkSpace
    sdg::SequenceDistanceGraph{LongDNASeq}
    paired_reads_datastores::Vector{PairedReads}
    long_reads_datastores::Vector{LongReads}
    linked_reads_datastores::Vector{LinkedReads}
end

function WorkSpace()
    return WorkSpace(SequenceDistanceGraph{LongDNASeq}(),
                     Vector{PairedReads}(),
                     Vector{LongReads}(),
                     Vector{LinkedReads}())
end

function add_paired_reads!(ws::WorkSpace, reads::PairedReads)
    push!(ws.paired_reads_datastores, reads)
    return ws
end

function add_paired_reads!(ws::WorkSpace, file::String)
    reads = open(PairedReads, file)
    add_paired_reads!(ws, reads)
    return ws
end

function Base.show(io::IO, ws::WorkSpace)
    println("Graph Genome workspace")
    
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