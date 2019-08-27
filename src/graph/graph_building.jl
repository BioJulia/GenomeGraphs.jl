

#=
function new_graph_from_kmerlist(kmerlist::Vector{DNAMer{K}}) where {K}
    str = string("onstructing Sequence Distance Graph from ", length(kmerlist), ' ', K, "-mers")
    @info string('C', str)
    sg = GRAPH_TYPE()
    build_unitigs_from_sorted_kmers!(sg, kmerlist)
    if n_nodes(sg) > 1
        connect_unitigs_by_overlaps!(sg, DNAMer{K})
    end
    @info string("Done c", str)
    return sg
end
SequenceDistanceGraph(kmerlist::Vector{DNAMer{K}}) where {K} = new_graph_from_kmerlist(kmerlist)
=#
