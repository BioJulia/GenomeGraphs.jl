struct SequenceGraphPath
    sg::SequenceDistanceGraph
    nodes::Vector{NodeID}
end

@inline nodes(p::SequenceGraphPath) = p.nodes
@inline push!(p::SequenceGraphPath, n::NodeID) = push!(nodes(p), n)

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