

struct SequenceGraph
    nodes::Vector{SequenceGraphNode}
    links::Vector{Vector{SequenceGraphLink}}
end

nodes(sg::SequenceGraph) = sg.nodes
node(sg::SequenceGraph, i::NodeID) = nodes(sg)[abs(i)] 
links(sg::SequenceGraph) = sg.links

"""
    links(sg::SequenceGraph, node::NodeID)

Get all of the links of a Node of a sequence graph.
"""
function links(sg::SequenceGraph, node::NodeID)
    l = links(sg)
    @inbounds return l[abs(node)]
end


"""
    add_node!(sg::SequenceGraph, n::SequenceGraphNode)

Add a SequenceGraphNode `n` to a SequenceGraph `sg`.

Returns the node ID used to access the node added from the graph.
"""
function add_node!(sg::SequenceGraph, n::SequenceGraphNode)
    newlen = length(push!(nodes(sg)))
    resize!(links(sg), newlen)
    return newlen
end

function remove_node!(sg::SequenceGraph, n::NodeID)
    oldlinks = copy(links(sg, n))
    for oldlink in oldlinks
        remove_link!(sg, source(oldlink), destination(oldlink))
    end
end


"""
    add_link!(sg::SequenceGraph, source::NodeID, dest::NodeID, dist::Int)

Construct a link between two nodes in a sequence Graph.
"""
function add_link!(sg::SequenceGraph, source::NodeID, dest::NodeID, dist::Int)
    push!(links(sg, source), Link(source, dest, dist))
    push!(links(sg, dest), Link(dest, source, dist))
end

function remove_link!(sg::SequenceGraph, source::NodeID, dest::NodeID)
    slinks = links(sg, source)
    slinkslen = length(slinks)
    filter!(!isequal(SequenceGraphLink(source, dest, 0)), slinks)
    dlinks = links(sg, dest)
    dlinkslen = length(dlinks)
    filter!(!isequal(SequenceGraphLink(dest, source, 0)), dlinks)
    return slinkslen != length(slinks) || dlinkslen != length(dlinks)
end






"""
    forward_links(sg::SequenceGraph, n::NodeID)

Get the list of links leaving `n`.
"""
function forward_links(sg::SequenceGraph, n::NodeID)
    r = Vector{Link}(0)
    nodelinks = links(sg, n)
    sizehint!(r, length(nodelinks))
    for link in nodelinks
        if is_forward_link(link, n)
            push!(r, link)
        end
    end
    return r
end

backward_links(sg::SequenceGraph, n::NodeID) = forward_links(sg, -n)
