
struct SequenceGraph
    nodes::Vector{Node}
    links::Vector{Vector{Link}}
end

nodes(sg::SequenceGraph) = sg.nodes
links(sg::SequenceGraph) = sg.links

"""
    aadd_node!(sg::SequenceGraph, n::SequenceGraphNode)

Add a SequenceGraphNode `n` to a SequenceGraph `sg`.

Returns the node ID used to access the node added from the graph.
"""
function add_node!(sg::SequenceGraph, n::SequenceGraphNode)
    push!(nodes(sg))
    linksref = links(sg)
    newlen = length(linksref) + 1
    resize!(linksref, newlen)
    return newlen
end

function remove_node!(sg::SequenceGraph, n::NodeID)
    oldlinks = copy(links(sg))

end


"""
    add_link!(sg::SequenceGraph, source::NodeID, dest::NodeID, dist::Int)

Construct a link between two nodes in a sequence Graph.
"""
function add_link!(sg::SequenceGraph, source::NodeID, dest::NodeID, dist::Int)
    ls = links(sg)
    push!(ls[abs(source)], Link(source, dest, dist))
    push!(ls[abs(dest)], Link(dest, source, dist))
end

void SequenceGraph::add_link(sgNodeID_t source, sgNodeID_t dest, int32_t d) {
    Link l(source,dest,d);
    links[(source > 0 ? source : -source)].emplace_back(l);
    std::swap(l.source,l.dest);
    links[(dest > 0 ? dest : -dest)].emplace_back(l);
}

function remove_link!(sg::SequenceGraph, )

end
