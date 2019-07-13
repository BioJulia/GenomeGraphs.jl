###
### New unitig graph construction straight from a set of kmers
###

# An improvement on the old algorithm:
# Builds a graph of unitigs straight away.
# Uses a sorted list (vector) of canonical kmers instead of a set. Checking kmer
# membership in the list requires binary lookup, which is slower than the O(1) of
# a set's hash, but typical sets of kmers for bigger genomes use waaay to much
# memory, so a list it is. But, this algorithm tries to query the list as few times
# as possible, and keeps track of all kmers used up in unitigs.
# Also, unlike the previous algorithm, when it came to making connections, an all
# v all double for loop was used, which was a bit pointless, so we avoid that
# in this version, since the two connection vectors are sorted, we can do a single
# pass of the two vectors instead. 
# With decent Kmer counting, this thing should be able to at least put together
# medium size genome like arabidopsis without much trouble. Larger genomes are
# probably fine too, but might need a big machine.
# If we use a large K value that will be allowed by the Skipmer types comming
# in BioSequences v2, the unitig graphs this will produce for medium size genomes
# should be pretty darned good in the first place.

struct Kidx{K}
    kmer::DNAKmer{K}
    idx::UInt64
end

encoded_data(mer) = reinterpret(UInt64, mer)

function iscanonical(seq)
    i = 1
    j = length(seq)
    @inbounds while i <= j
        f = seq[i]
        r = complement(seq[j])
        f < r && return true
        r < f && return false
        i += 1
        j -= 1
    end
    return true
end

function canonical!(seq::BioSequence{DNAAlphabet{2}})
    if !iscanonical(seq)
        reverse_complement!(seq)
    end
    return seq
end

# TODO: Update BioSequences with new neighbour iterators instead. For now, these
# functions will do.
function kmer_fw_neighbours(mer::DNAKmer{K}) where {K}
    d = encoded_data(mer)
    base = d << 2
    return (DNAKmer{K}(base), DNAKmer{K}(base + 0x01), DNAKmer{K}(base + 0x02), DNAKmer{K}(base + 0x03))
end

function kmer_bw_neighbours(mer::DNAKmer{K}) where {K}
    d = encoded_data(mer)
    base = d >> 2
    BT = typeof(base)
    return (DNAKmer{K}(base), DNAKmer{K}(base + (BT(1) << 2(K - 1))), DNAKmer{K}(base + (BT(2) << 2(K - 1))), DNAKmer{K}(base + (BT(3) << 2(K - 1))))
end

function is_end_bw(mer::DNAKmer{K}, merlist::Vector{DNAKmer{K}}) where {K}
    @debug "Checking if kmer is end BW" mer
    next = Vector{Kidx{K}}()
    get_bw_idxs!(next, mer, merlist)
    @debug "BW neighbours:" next
    length(next) != 1 && return true
    @inbounds p = next[1].kmer
    get_fw_idxs!(next, p, merlist)
    @debug "FW neighbours of only BW neighbour:" p next
    length(next) != 1 && return true
    return false
end

function is_end_fw(mer::DNAKmer{K}, merlist::Vector{DNAKmer{K}}) where {K}
    @debug "Checking if kmer is end FW" mer
    next = Vector{Kidx{K}}()
    get_fw_idxs!(next, mer, merlist)
    @debug "FW neighbours:" next
    length(next) != 1 && return true
    @inbounds p = next[1].kmer
    get_bw_idxs!(next, p, merlist)
    @debug "BW neighbours of only FW neighbour:" p next
    length(next) != 1 && return true
    return false
end

function get_bw_idxs!(out::Vector{Kidx{K}}, kmer::DNAKmer{K}, kmerlist::Vector{DNAKmer{K}}) where {K}
    empty!(out)
    for n in kmer_bw_neighbours(kmer)
        cnext = canonical(n)
        cidx = min(searchsortedfirst(kmerlist, cnext), length(kmerlist))
        if @inbounds kmerlist[cidx] == cnext
            push!(out, Kidx{K}(n, cidx))
        end
    end
end

function get_fw_idxs!(out::Vector{Kidx{K}}, kmer::DNAKmer{K}, kmerlist::Vector{DNAKmer{K}}) where {K}
    empty!(out)
    for n in kmer_fw_neighbours(kmer)
        cnext = canonical(n)
        cidx = min(searchsortedfirst(kmerlist, cnext), length(kmerlist))
        if @inbounds kmerlist[cidx] == cnext
            push!(out, Kidx{K}(n, cidx))
        end
    end
end

const GRAPH_TYPE = SequenceDistanceGraph{BioSequence{DNAAlphabet{2}}}

function build_unitigs_from_kmerlist!(sg::GRAPH_TYPE, kmerlist::Vector{DNAKmer{K}}) where {K}
    @info string("Constructing unitigs from ", length(kmerlist), " ", K, "-mers")
    used_kmers = falses(length(kmerlist))
    
    for start_kmer_idx in eachindex(kmerlist)
        @debug "Considering new kmer" start_kmer_idx
        
        # Any kmer can only occur in one unitig.
        if used_kmers[start_kmer_idx]
            @debug "Kmer has been used" start_kmer_idx
            continue
        end
        
        # Check if the kmer is an end/junction of a unitig.
        start_kmer = kmerlist[start_kmer_idx]
        end_bw = is_end_bw(start_kmer, kmerlist)
        end_fw = is_end_fw(start_kmer, kmerlist)
        
        if !end_bw && !end_fw
            @debug "Kmer is middle of a unitig" start_kmer_idx start_kmer
            continue
        end
        
        if end_bw && end_fw
            @debug "Kmer is single unitig" start_kmer_idx start_kmer 
            # Kmer as unitig
            s = BioSequence{DNAAlphabet{2}}(start_kmer)
            used_kmers[start_kmer_idx] = true
        else
            # A unitig starts on this kmer.
            # Make sure the unitig starts on FW.
            current_kmer = start_kmer
            used_kmers[start_kmer_idx] = true
            if end_fw
                current_kmer = reverse_complement(start_kmer)
                end_fw = end_bw
            end
            @debug "Start of unitig" start_kmer current_kmer end_bw end_fw
            # Start unitig 
            s = BioSequence{DNAAlphabet{2}}(current_kmer)
            fwn = Vector{Kidx{K}}()
            while !end_fw
                # Add end nucleotide, update current kmer.
                get_fw_idxs!(fwn, current_kmer, kmerlist)
                @debug "Extending unitig" fwn 
                current_kmer = first(fwn).kmer
                if used_kmers[first(fwn).idx]
                    @debug "New kmer is already used" current_kmer
                    break # Break circular contigs into lines.
                end
                used_kmers[first(fwn).idx] = true
                push!(s, last(current_kmer))
                end_fw = is_end_fw(current_kmer, kmerlist)
            end
        end
        add_node!(sg, canonical!(s))
    end
    # A temporary check for circle problem for now.
    if !all(used_kmers)
        @warn "Some kmers have not been incorporated into unitigs. This may be a case of the circle problem" kmerlist[(!).(used_kmers)]
    end
    @info string("Constructed ", length(nodes(sg)), " unitigs")
    return sg
end

function find_unitig_overlaps(sg::GRAPH_TYPE, ::Type{DNAKmer{K}}) where {K}
    @info string("Identifying the ", K - 1, "bp (K - 1) overlaps between ", length(nodes(sg)), " unitigs")
    # Save the (k-1)mer in (rev on first k-1 / fw on last k-1) or out ( fw on first k-1 / bw on last k-1)
    @debug "Sorting K - 1 overlaps as `in` or `out`"
    in = Vector{Tuple{DNAKmer{K-1},NodeID}}()
    out = Vector{Tuple{DNAKmer{K-1},NodeID}}()
    sizehint!(in, length(nodes(sg)))
    sizehint!(out, length(nodes(sg)))
    for nid in eachindex(nodes(sg))
        nodeseq = node(sg, nid).seq
        firstmer = DNAKmer{K-1}(nodeseq[1:K - 1])
        @debug string("Considering node ", nid) nodeseq
        if iscanonical(firstmer)
            @debug "Source overlap is canonical"
            push!(in, (firstmer, nid))
        else
            @debug "Source overlap is not canonical"
            push!(out, (reverse_complement(firstmer), nid))
        end
        lastmer = DNAKmer{K-1}(nodeseq[end - (K - 2):end])
        if iscanonical(lastmer)
            @debug "Sink overlap is canonical"
            push!(out, (lastmer, -nid))
        else
            @debug "Sink overlap is not canonical"
            push!(in, (reverse_complement(lastmer), -nid))
        end
    end
    sort!(in)
    sort!(out)
    return in, out
end

function connect_unitigs_by_overlaps!(sg::GRAPH_TYPE, ::Type{DNAKmer{K}}) where {K}
    in, out = find_unitig_overlaps(sg, DNAKmer{K})
    ol = length(out)
    @info string("Linking ", length(nodes(sg)), " unitigs by their ", K - 1, "bp (K - 1) overlaps")
    # Connect all out -> in for all combinations on each kmer.
    next_out_idx = 1
    for i in in
        while next_out_idx <= ol && first(out[next_out_idx]) < first(i)
            next_out_idx += 1
        end
        oidx = next_out_idx
        while oidx <= ol && first(out[oidx]) == first(i)
            add_link!(sg, last(i), last(out[oidx]), -K + 1) # No support, although we could add the DBG operation as such.
            oidx += 1
        end
    end
end

function new_graph_from_kmerlist(kmerlist::Vector{DNAKmer{K}}) where {K}
    str = string("onstructing Sequence Distance Graph from ", length(kmerlist), ' ', K, "-mers")
    @info string('C', str)
    sg = GRAPH_TYPE()
    build_unitigs_from_kmerlist!(sg, kmerlist)
    if n_nodes(sg) > 1
        connect_unitigs_by_overlaps!(sg, DNAKmer{K})
    end
    @info string("Done c", str)
    return sg
end
SequenceDistanceGraph(kmerlist::Vector{DNAKmer{K}}) where {K} = new_graph_from_kmerlist(kmerlist)

