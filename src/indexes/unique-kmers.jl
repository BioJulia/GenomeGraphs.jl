struct GraphStrandPosition
    node::NodeID
    position::UInt32
end

Base.summary(io::IO, pos::GraphStrandPosition) = print(io, "(Node: ", pos.node, ", Base position: ", pos.position, ')')
Base.show(io::IO, pos::GraphStrandPosition) = summary(io, pos)

struct UniqueMerIndex{M<:AbstractMer}
    mer_to_graphposition::Dict{M,GraphStrandPosition}
    unique_mers_per_node::Vector{UInt64}
    total_mers_per_node::Vector{UInt64}
end

@inline function isunmappable(idx::UniqueMerIndex, nodeid::NodeID)
    return idx.unique_mers_per_node[abs(nodeid)] == 0
end

function find_unique_mer(idx::UniqueMerIndex{M}, mer::M) where {M<:AbstractMer}
    pos = get(idx.mer_to_graphposition, mer, GraphStrandPosition(0, 0))
    exists = pos != GraphStrandPosition(0, 0)
    return exists, pos
end

function _discard_nonunique_mers!(mers::Vector{Pair{M,GraphStrandPosition}}) where {M<:AbstractMer}
    wi = 1
    ri = 1
    nri = 1
    stop = lastindex(mers) + 1
    while ri < stop
        while nri != stop && first(mers[nri]) == first(mers[ri])
            nri = nri + 1
        end
        if nri - ri == 1
            mers[wi] = mers[ri]
            wi = wi + 1
        end
        ri = nri
    end
    resize!(mers, wi - 1)
    return mers
end

function UniqueMerIndex{M}(graph::S) where {S<:SequenceDistanceGraph,M<:AbstractMer}
    t = 0
    total_mers_per_node = zeros(UInt64, n_nodes(graph))
    for nodeid in each_node_id(graph)
        seq = sequence(graph, nodeid)
        if length(seq) >= BioSequences.ksize(M)
            n = length(seq) + 1 - BioSequences.ksize(M)
            t = t + n
            total_mers_per_node[nodeid] = n
        end
    end 
    kidxv = Vector{Pair{M,GraphStrandPosition}}(undef, t)
    index_i = 1
    for nodeid in each_node_id(graph)
        seq = sequence(graph, nodeid)
        if length(seq) >= BioSequences.ksize(M)
            for mer in each(M, seq)
                merpos = BioSequences.position(mer)
                mernode = ifelse(fwmer(mer) <= bwmer(mer), nodeid, -nodeid)
                mermer = ifelse(fwmer(mer) <= bwmer(mer), fwmer(mer), bwmer(mer))
                kidxv[index_i] = Pair(mermer, GraphStrandPosition(mernode, merpos))
                index_i = index_i + 1
            end
        end
    end
    
    sort!(kidxv, by = x -> first(x))
    
    _discard_nonunique_mers!(kidxv)
    
    mer_to_graphposition = Dict{M, GraphStrandPosition}()
    sizehint!(mer_to_graphposition, length(kidxv))
    unique_mers_per_node = zeros(UInt64, length(kidxv))
    
    for kidx in kidxv
        mer_to_graphposition[first(kidx)] = last(kidx)
        unique_mers_per_node[abs(last(kidx).node)] += 1
    end
    
    return UniqueMerIndex{M}(mer_to_graphposition, unique_mers_per_node, total_mers_per_node)
end
