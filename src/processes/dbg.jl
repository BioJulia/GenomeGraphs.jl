###
### de-bruijn graph and initial workspace construction
###

# An improvement on the old algorithm:
# Builds a graph of unitigs straight away.
# Uses a sorted list (vector) of canonical kmers instead of a set. Checking kmer
# membership in the list requires binary lookup, which is slower than the O(1) of
# a set's hash, but typical sets of kmers for bigger genomes use waaaaay to much
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

struct Kidx{M<:AbstractMer}
    kmer::M
    idx::UInt64
end

# TODO: Update BioSequences with new neighbour iterators instead. For now, these
# functions will do.
function kmer_fw_neighbours(mer::M) where {M<:AbstractMer}
    d = BioSequences.encoded_data(mer)
    base = d << 2
    return (M(base), M(base + 0x01), M(base + 0x02), M(base + 0x03))
end

function kmer_bw_neighbours(mer::M) where {M<:AbstractMer}
    k = BioSequences.ksize(M) - 1
    d = BioSequences.encoded_data(mer)
    base = d >> 2
    BT = typeof(base)
    return (M(base), M(base + (BT(1) << 2(k))), M(base + (BT(2) << 2(k))), M(base + (BT(3) << 2(k))))
end

function is_end_bw(mer::M, merlist::Vector{M}) where {M<:AbstractMer}
    @debug "Checking if kmer is end BW" mer
    next = Vector{Kidx{M}}()
    get_bw_idxs!(next, mer, merlist)
    @debug "BW neighbours:" next
    length(next) != 1 && return true
    @inbounds p = next[1].kmer
    get_fw_idxs!(next, p, merlist)
    @debug "FW neighbours of only BW neighbour:" p next
    length(next) != 1 && return true
    return false
end

function is_end_fw(mer::M, merlist::Vector{M}) where {M<:AbstractMer}
    @debug "Checking if kmer is end FW" mer
    next = Vector{Kidx{M}}()
    get_fw_idxs!(next, mer, merlist)
    @debug "FW neighbours:" next
    length(next) != 1 && return true
    @inbounds p = next[1].kmer
    get_bw_idxs!(next, p, merlist)
    @debug "BW neighbours of only FW neighbour:" p next
    length(next) != 1 && return true
    return false
end

function get_bw_idxs!(out::Vector{Kidx{M}}, kmer::M, kmerlist::Vector{M}) where {M<:AbstractMer}
    empty!(out)
    for n in kmer_bw_neighbours(kmer)
        cnext = canonical(n)
        cidx = min(searchsortedfirst(kmerlist, cnext), length(kmerlist))
        if @inbounds kmerlist[cidx] == cnext
            push!(out, Kidx{M}(n, cidx))
        end
    end
end

function get_fw_idxs!(out::Vector{Kidx{M}}, kmer::M, kmerlist::Vector{M}) where {M<:AbstractMer}
    empty!(out)
    for n in kmer_fw_neighbours(kmer)
        cnext = canonical(n)
        cidx = min(searchsortedfirst(kmerlist, cnext), length(kmerlist))
        if @inbounds kmerlist[cidx] == cnext
            push!(out, Kidx{M}(n, cidx))
        end
    end
end

const GRAPH_TYPE = SequenceDistanceGraph{LongSequence{DNAAlphabet{4}}}

function build_unitigs_from_sorted_kmers!(
    sg::SequenceDistanceGraph{LongSequence{A}},
    kmerlist::Vector{M}) where {A<:DNAAlphabet,M<:AbstractMer{DNAAlphabet{2}}}
    
    
    @info string("Constructing unitigs from ", length(kmerlist), " ", BioSequences.ksize(M), "-mers")
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
            s = LongSequence{A}(start_kmer)
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
            s = LongSequence{A}(current_kmer)
            fwn = Vector{Kidx{M}}()
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

minus_one_k(::Type{Mer{A,K}}) where {A,K} = Mer{A,K-1}
minus_one_k(::Type{BigMer{A,K}}) where {A,K} = BigMer{A,K-1}

function find_unitig_overlaps(sg::GRAPH_TYPE, ::Type{M}) where {M<:AbstractMer}
    @info string("Identifying the ", BioSequences.ksize(M) - 1, "bp (K - 1) overlaps between ", length(nodes(sg)), " unitigs")
    # Save the (k-1)mer in (rev on first k-1 / fw on last k-1) or out ( fw on first k-1 / bw on last k-1)
    @debug "Sorting K - 1 overlaps as `in` or `out`"
    
    in = Vector{Tuple{minus_one_k(M),NodeID}}()
    out = Vector{Tuple{minus_one_k(M),NodeID}}()
    sizehint!(in, length(nodes(sg)))
    sizehint!(out, length(nodes(sg)))
    for nid in eachindex(nodes(sg))
        nodeseq = node(sg, nid).seq
        firstmer = minus_one_k(M)(nodeseq[1:BioSequences.ksize(M) - 1])
        @debug string("Considering node ", nid) nodeseq
        if iscanonical(firstmer)
            @debug "Source overlap is canonical"
            push!(in, (firstmer, nid))
        else
            @debug "Source overlap is not canonical"
            push!(out, (reverse_complement(firstmer), nid))
        end
        lastmer = minus_one_k(M)(nodeseq[end - (BioSequences.ksize(M) - 2):end])
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

function connect_unitigs_by_overlaps!(sg::GRAPH_TYPE, ::Type{M}) where {M<:AbstractMer}
    in, out = find_unitig_overlaps(sg, M)
    ol = length(out)
    @info string("Linking ", length(nodes(sg)), " unitigs by their ", BioSequences.ksize(M) - 1, "bp (K - 1) overlaps")
    # Connect all out -> in for all combinations on each kmer.
    next_out_idx = 1
    for i in in
        while next_out_idx <= ol && first(out[next_out_idx]) < first(i)
            next_out_idx += 1
        end
        oidx = next_out_idx
        while oidx <= ol && first(out[oidx]) == first(i)
            add_link!(sg, last(i), last(out[oidx]), -BioSequences.ksize(M) + 1) # No support, although we could add the DBG operation as such.
            oidx += 1
        end
    end
end

"""
"""
function dbg(::Type{M}, min_freq::Integer, file::String, name::Union{String,Nothing} = nothing) where {M<:AbstractMer}
    ws = WorkSpace()
    ds = open(PairedReads, file, name)
    add_paired_reads!(ws, ds)
    _dbg!(ws.sdg, ds, M, UInt8(min_freq))
    return ws
end

"""
    dbg!(ws::WorkSpace, ds::String, ::Type{M}, min_freq::Integer) where {M<:AbstractMer}


"""
function dbg!(ws::WorkSpace, ::Type{M}, min_freq::Integer, name::String) where {M<:AbstractMer}
    reads = paired_reads(ws, name)
    _dbg!(ws.sdg, reads, M, UInt8(min_freq))
    return ws
end

function _dbg!(sg::GRAPH_TYPE, ds::PairedReads, ::Type{M}, min_freq::UInt8) where {M<:AbstractMer}
    @info "Counting kmers in datastore"
    # In the future do a better kmer counting - but this will do for e.coli to prove a point.
    spectra = build_freq_list(M, buffer(ds), 1:Int(length(ds)))
    filter!(x -> freq(x) ≥ min_freq, spectra)
    merlist = [mer(x) for x in spectra]
    return _dbg!(sg, merlist)
end

function _dbg!(sg::GRAPH_TYPE, kmerlist::Vector{M}) where {M<:AbstractMer}
    str = string("onstructing compressed de-bruijn graph from ", length(kmerlist), ' ', BioSequences.ksize(M), "-mers")
    @info string('C', str)
    build_unitigs_from_sorted_kmers!(sg, kmerlist)
    if n_nodes(sg) > 1
        connect_unitigs_by_overlaps!(sg, M)
    end
    @info string("Done c", str)
    return sg
end