struct SequenceGraphPath{G<:SequenceDistanceGraph}
    sg::G
    nodes::Vector{NodeID}
end

SequenceGraphPath(sg::G) where {G<:SequenceDistanceGraph} = SequenceGraphPath{G}(sg, Vector{NodeID}())

@inline nodes(p::SequenceGraphPath) = p.nodes
@inline n_nodes(p::SequenceGraphPath) = length(nodes(p))
@inline Base.push!(p::SequenceGraphPath, n::NodeID) = push!(nodes(p), n)
@inline graph(p::SequenceGraphPath) = p.sg
@inline Base.first(p::SequenceGraphPath) = first(nodes(p))
@inline Base.last(p::SequenceGraphPath) = last(nodes(p))

#=
function Base.reverse!(p::SequenceGraphPath)
    nds = nodes(p)
    for i in 1:div(lastindex(nds), 2)
        j = lastindex(nds) - i + 1
        @inbounds x = nds[i]
        @inbounds y = nds[j]
        @inbounds nds[i] = -y
        @inbounds nds[j] = -x
    end
    return p
end
=#

function Base.reverse!(p::SequenceGraphPath)
    nds = nodes(p)
    i = firstindex(nds)
    j = lastindex(nds)
    while i ≤ j
        @inbounds x = nds[i]
        @inbounds y = nds[j]
        @inbounds nds[i] = -y
        @inbounds nds[j] = -x
        i = i + 1
        j = j - 1
    end
    return p
end

# TODO: See about potential performance improvement using BioSequences.jl `mismatches`.
function check_overlap(a::BioSequence, b::BioSequence, by::Integer)
    i′ = lastindex(a) + 1
    i = i′ - by
    j = 1
    while i < i′ 
        if a[i] != b[j]
            return false
        end
        i = i + 1
        j = j + 1
    end
    return true
end

# TODO: This function currently does not check for deleted nodes in a graph.
# It probably should throw an error if a node is deleted and you're trying to
# use it in a path.
function sequence(p::SequenceGraphPath{SequenceDistanceGraph{S}}) where {S<:BioSequence}
    s = S()
    pnode = 0
    for n in nodes(p)
        if n > 0
            nseq = sequence(graph(p), n)
        else
            nseq = reverse_complement(sequence(graph(p), n))
        end
        if pnode != 0
            # Find link between pnode output (+pnode) and n's sink (-n)
            lnk = find_link(graph(p), pnode, n)
            if !isnothing(lnk)
                dst = distance(lnk)
                if dst > 0
                    
                else
                    # Check that the overlap is valid in terms of it's sequence.
                    ovl = -dst
                    if !check_overlap(s, nseq, ovl)
                        error("Sequences do not overlap.")
                    end
                    # TODO: Investigate if this can be done more efficiently
                    # using resize! and copyto!.... Or perhaps something more
                    # hacky.
                    nseq = nseq[ovl + 1:end]
                end
            else
                error("This path is not a valid path.")
            end
        end
        #oldlen = length(s)
        #resize!(s, oldlen + length(nseq) - ovl)
        #copyto!(s, oldlen + 1, nseq, ovl + 1, length(nseq) - ovl)
        append!(s, nseq)
        pnode = -n
    end
    return s
end

function join_path!(p::SequenceGraphPath, consume::Bool)
    pnodes = Set{NodeID}()
    for n in nodes(p)
        push!(pnodes, n)
        push!(pnodes, -n)
    end
    # If p is not canonical, flip it.
    ps = sequence(p)
    if !iscanonical(ps)
        reverse_complement!(ps)
        reverse!(p)
    end
    # Add a new node that is this sequence to the graph.
    newid = add_node!(graph(p), ps)
    # Connect the new node up.
    # TODO: Circular structures in the graph MIGHT screw this up, but I'm not sure.
    for l in backward_links(graph(p), first(p))
        add_link!(graph(p), newid, destination(l), distance(l))
    end
    for l in forward_links(graph(p), last(p))
        add_link!(graph(p), -newid, destination(l), distance(l))
    end
    
    if consume
        for n in pnodes
            # Check if the node has neighbours that are not involved in the path.
            ext_neigh = false
            if n != last(p)
                for l in forward_links(graph(p), n)
                    if destination(l) ∉ pnodes
                        ext_neigh = true
                    end
                end
            end
            if n != first(p)
                for l in backward_links(graph(p), n)
                    if destination(l) ∉ pnodes
                        ext_neigh = true
                    end
                end
            end
            if ext_neigh
                continue
            end
            remove_node!(graph(p), n)
        end
    end
    return newid
end