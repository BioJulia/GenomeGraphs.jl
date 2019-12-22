
function remove_tips!(ws::WorkSpace, min_size::Integer)
    @info "Beginning tip removal process"
    sdg = graph(ws)
    pass = 1
    tips = Graphs.find_tip_nodes(sdg, min_size)
    ntips = length(tips)
    utgs = Vector{Graphs.SequenceGraphPath{typeof(sdg)}}()
    newnodes = Vector{Graphs.NodeID}()
    while true
        @info string("Pass number: ", pass)
        if isempty(tips)
            @info "No more tips in graph"
            break
        end
        @info string("Found ", ntips, " tips")
        for tip in tips
            Graphs.remove_node!(sdg, tip)
        end
        @info string("Removed ", ntips, " tips")
        @info "Collapsing any resulting transient paths"
        Graphs.collapse_all_unitigs!(utgs, newnodes, sdg, 2, true)
        pass = pass + 1
        tips = Graphs.find_tip_nodes(sdg, min_size)
        ntips = length(tips)
    end
    @info string("Finished tip removal process in ", pass, " passes")
    return ws
end

"""
For sequences longer than initial kmers we can not naively check the overlaps between them
using the same functions designed for Kmers. This is because direction becomes an issue.

In order to be a backward neighbor the overlap must be between the reverse complement of the prefix of a sequence
or the suffix of the sequence itself. Otherwise the overlap is contained inside the sequence and is not a valid overlap.

Thus for each edge we produce 2 candidates that can be considered for backward neighbor search

"""
function backward_candidates(sdg::Graphs.SequenceDistanceGraph,K::Int64)
    cands = Vector{Tuple{DNAKmer{K},Int64}}()
    l = Base.length(nodes(sdg))
    nodes_ = nodes(sdg)
    for i in eachindex(nodes(sdg))
        if nodes(sdg)[i].deleted==true
            continue
        end
        kpref = DNAKmer{K}(reverse_complement(sequence(nodes_[i])[1:K]))
        ksuff = DNAKmer{K}(sequence(nodes_[i])[Base.length(sequence(nodes_[i]))-K+1:end])
        push!(cands,(kpref,-i))
        push!(cands,(ksuff,i))
    end
    return cands
end

function forward_candidates(sdg::Graphs.SequenceDistanceGraph,K::Int64)
    cands = Vector{Tuple{DNAKmer{K},Int64}}()
    l = Base.length(nodes(sdg))
    nodes_ = nodes(sdg)
    for i in eachindex(nodes(sdg))
        if nodes(sdg)[i].deleted==true
            continue
        end
        kpref = DNAKmer{K}(sequence(nodes_[i])[1:K])
        ksuff = DNAKmer{K}(reverse_complement(sequence(nodes_[i])[Base.length(sequence(nodes_[i]))-K+1:end]))
        push!(cands,(kpref,i))
        push!(cands,(ksuff,-i))
    end
    return cands
end

function delete_tips_graph(sdg::Graphs.SequenceDistanceGraph,coverage :: Vector{Float64},K::Int64)
    @info string("Coverage length : " , Base.length(coverage))
    @info string("Node length : " , Base.length(nodes(sdg)))
    @assert Base.length(nodes(sdg)) == Base.length(coverage)  "Check coverage information"
    nodes_ = nodes(sdg)
    tip_deleted = false
    used_seqs = falses(Base.length(nodes_))
    prefs = get_prefixes(sdg,K)
    suffs = get_suffixes(sdg,K)
    forward_cands_with_ids   = forward_candidates(sdg,K)
    backward_cands_with_ids  = backward_candidates(sdg,K)
    sort!(forward_cands_with_ids)
    sort!(backward_cands_with_ids)
    #back_cands = [x[1] for x in backward_cands_with_ids]
    #forw_cands = [x[1] for x in forward_cands_with_ids]
    parents = Vector{Int64}()
    tips = Vector{Tuple{Int64,Int64}}()
    all_kmers_with_ids = vcat(prefs,suffs)
    sort!(all_kmers_with_ids)
    all_kmers = [tup[1] for tup in all_kmers_with_ids] ## get kmers from (kmer,index) tuples
    #@info string("All $K-mers: ", all_kmers)
    @info string("Forward candidates " , forward_cands_with_ids)
    @info string("Backward candidates " , backward_cands_with_ids)
    for start_seq_idx in eachindex(nodes_)
        if nodes_[start_seq_idx].deleted
            @info string("Skipping deleted node ", start_seq_idx)
            continue
        end
        #start_sequence_pref = prefs[start_seq_idx]
        #start_sequence_suff = suffs[start_seq_idx]
        pref_seq = get_prefix(nodes_[start_seq_idx].seq,K)
        rev_pref_seq = reverse_complement(pref_seq)
        suff_seq = get_suffix(nodes_[start_seq_idx].seq,K)
        rev_suff_seq = reverse_complement(suff_seq)
        next = Vector{Kidx{K}}()
        next2 = Vector{Kidx{K}}()
        get_bw_idxs2!(next,pref_seq,backward_cands_with_ids)
        get_fw_idxs2!(next,rev_pref_seq,forward_cands_with_ids)
        get_fw_idxs2!(next2,suff_seq,forward_cands_with_ids)
        get_bw_idxs2!(next2,rev_suff_seq,backward_cands_with_ids)
        #next2 = Vector{Kidx{K}}() # not used
        #get_bw_idxs!(next2, start_sequence_pref, all_kmers)
        @info string("Backward neighbors of ", pref_seq, " is: ", next)
        @info string("Forward  neighbors of ", suff_seq, " is: ", next2)
        if (Base.length(next)==0 && Base.length(next2)==1) || (Base.length(next)==1 && Base.length(next2)==0)
            pars_ = vcat(next,next2)
            parent = first(pars_).idx
            @info string(nodes_[start_seq_idx] , " is a tip contig!! with parent: ", parent)
            push!(tips,(parent,start_seq_idx))
            push!(parents,parent)
        end
    end
    sort!(tips)
    sort!(parents)
    i = 1
    @info string("All found tips : ",tips)
    @info string("Parents : ", parents)
    while i <= Base.length(tips)
        clust = searchsorted(parents,tips[i][1])
        @info string("Tip number : ", tips[i])
        @info string("Clust :  ", clust)
        if Base.length(clust)==1
            i = clust[end]+1
            continue
            start_ind = clust[1]
            end_ind   = clust[1]
        else
            start_ind = clust[1]
            end_ind   = clust[2]
        end
        real_inds = [tips[x][2] for x in clust ]
        max_ind  = real_inds[argmax(map(Base.length,coverage[real_inds]))]
        for ind in real_inds
            if ind!=max_ind
                remove_node!(sdg,ind)
                tip_deleted = true
            end
        end
        i  = clust[end]+1
    end
    @info string("After tip removal, ", sdg)
    return sdg, tip_deleted
end