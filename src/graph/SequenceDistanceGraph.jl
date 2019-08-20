###
### Node and Link types for SequenceDistanceGraph
###

const NodeID = Int64

include("SDGNode.jl")
include("DistanceGraphLink.jl")

const LinksT = Vector{Vector{DistanceGraphLink}}

###
### Graph types
###

"""
The SequenceDistanceGraph is a representation of a genome assembly.
Sequences are contained in nodes, and the distances are represented by links.

A singe node represents a sequence *and* its reverse complement.
Every node has a correlative ID starting from 1.
For every node X in the graph, the negative ID -X is mapped to the reverse
complement of X. This mapping is virtual: Only one node is stored in the graph.
This is because every node has an orientaton: Each node has a positive end (+),
and a negative end (-).
So when a node is accessed with (or traversed by entering) the positive end the
node yields the stored sequence.
Conversely, when a node is accessed with (or traversed by entering) the negative
end the node yelds the reverse complement of the stored sequence.
In this way the positive end can be thought of as the sequence start, and the
negative end can be thought of as the sequence end.

A single distance between two sequences is represented as a single link.
Every link connects two node ends and contains a distance (they take the form
`([+, -]n1, [+, -]n2, [+, -]dist)`).
A link connects two node ends, and so the order of the signed nodes in the links
does not change the link.
If the distance in a link is negative, this represents an overlap between two
sequences. These overlaps must be "perfect overlaps".
"""
struct SequenceDistanceGraph{S<:BioSequence}
    nodes::Vector{SDGNode{S}}
    links::LinksT
end

include("SequenceGraphPath.jl")

"Construct an empty sequence distance graph."
function SequenceDistanceGraph{S}() where {S<:BioSequence}
    return SequenceDistanceGraph{S}(Vector{SDGNode{S}}(), LinksT())
end

###
### Basic query and property access functions
###

"""
Get a reference to the vector of nodes in a graph `sg`.

!!! warning
    It is a bad idea to edit this vector yourself unless you know what you are
    doing. 
"""
@inline nodes(sg::SequenceDistanceGraph) = sg.nodes

"Get the number of nodes in the sequence distance graph `sg`."
@inline n_nodes(sg::SequenceDistanceGraph) = length(nodes(sg))

"Iterate over every node ID in the sequence distance graph `sg`."
@inline each_node_id(sg::SequenceDistanceGraph) = eachindex(nodes(sg))

"""
Get a reference to the vector of vectors of links in a graph `sg`.

!!! warning
    It is a bad idea to edit this vector yourself unless you know what you are doing.
"""
@inline links(sg::SequenceDistanceGraph) = sg.links

@inline function check_node_id(sg::SequenceDistanceGraph, i::NodeID)
    if 0 < abs(i) ≤ n_nodes(sg)
        return true
    end
    @error "Sequence graph has no node with ID of $i"
end

@inline node_unsafe(sg::SequenceDistanceGraph, n::NodeID) = @inbounds nodes(sg)[abs(n)]

"""
    node(sg::SequenceDistanceGraph, n::NodeID)

Get a specific node from a sequence distance graph `sg` using its
correlative node id `n`.

!!! note
    `node` accepts a NodeID that can be positive or negative.
    E.g. providing either 5 or -5 both mean node 5 in a graph,
    and so you will get the links for node 5.
"""
@inline function node(sg::SequenceDistanceGraph, n::NodeID)
    check_node_id(sg, n)
    return node_unsafe(sg, n)
end

@inline links_unsafe(sg::SequenceDistanceGraph, n::NodeID) = @inbounds links(sg)[abs(n)]

"""
    links(sg::SequenceGraph, n::NodeID)

Get all of the links of a Node of a sequence distance graph using its 
correlative node id `n`.

!!! note
    `links` accepts a NodeID that can be positive or negative.
    E.g. providing either 5 or -5 both mean node 5 in a graph,
    and so you will get the links for node 5.
"""
@inline function links(sg::SequenceDistanceGraph, n::NodeID)
    check_node_id(sg, n)
    return links_unsafe(sg, n)
end

@inline sequence_unsafe(sg::SequenceDistanceGraph, n::NodeID) = sequence(node_unsafe(sg, n))

"""
    sequence(sg::SequenceDistanceGraph, n::NodeID)

Get the full sequence of a node in a sequence distance graph using its
correlative node id `n`.
"""
function sequence(sg::SequenceDistanceGraph, n::NodeID)
    check_node_id(sg, n)
    return sequence_unsafe(sg, n)
end


###
### Graph editing operations
###

"""
    add_node!(sg::SequenceDistanceGraph{S}, n::SDGNode{S}) where {S<:BioSequence}

Add a node to a sequence distance graph.

Returns the node ID used to access the new node added in the graph.

!!! warning
    Currently, we don't enforce the sequence in the node is canonical here.
    We just trust that it is canonical.

!!! note
    Adding a node to the graph does just that. After adding the node it still
    will not be linked to any other nodes.
"""
function add_node!(sg::SequenceDistanceGraph{S}, n::SDGNode{S}) where {S<:BioSequence}
    newlen = length(push!(nodes(sg),n))
    push!(links(sg), Vector{DistanceGraphLink}())
    return newlen
end

"""
   add_node!(sg::SequenceDistanceGraph{S}, seq::BioSequence) where {S<:BioSequence}

Add a sequence to a sequence distance graph as a node.

Returns the node ID used to access the new node added in the graph.

Can accept any sequence type and will attempt to coerce the input sequence to the
type required by the graph.

!!! warning
    Currently, we don't enforce the sequence in the node is canonical here.
    We just trust that it is canonical.

!!! note
    Adding a node to the graph does just that. After adding the node it still
    will not be linked to any other nodes.
"""
function add_node!(sg::SequenceDistanceGraph{S}, seq::BioSequence) where {S<:BioSequence}
    return add_node!(sg, SDGNode{S}(convert(S, seq), false))
end

"""
    remove_node!(sg::SequenceDistanceGraph{S}, n::NodeID) where {S<:BioSequence}

Remove a node from a sequence distance graph.

!!! note
    Links involving this node will also be removed from the graph.
"""
function remove_node!(sg::SequenceDistanceGraph{S}, n::NodeID) where {S<:BioSequence}
    oldlinks = copy(links(sg, n))
    for oldlink in oldlinks
        remove_link!(sg, source(oldlink), destination(oldlink))
    end
    # TODO: This is a lazy solution to getting rid of the node.
    nodes[n] = empty_node(S)
end


"""
    add_link!(sg::SequenceDistanceGraph, source::NodeID, dest::NodeID, dist::Int)

Construct a link between two nodes in a sequence Graph.
"""
function add_link!(sg::SequenceDistanceGraph, source::NodeID, dest::NodeID, dist::Int)
    # Guard against someone adding links using ID's bigger than the current max NodeID in the graph.
    if abs(source) > length(links(sg))
        resize!(links(sg), abs(source))
    end
    if abs(dest) > length(links(sg))
        resize!(links(sg), abs(dest))
    end
    push!(links(sg, source), DistanceGraphLink(source, dest, dist))
    push!(links(sg, dest), DistanceGraphLink(dest, source, dist))
end

"""
    remove_link!(sg::SequenceDistanceGraph, source::NodeID, dest::NodeID)

Remove a link between two nodes in a SequenceDistanceGraph.
Returns a boolean indicating whether the removal was successful.
Reasons this function would not return `true` include that the link
didn't exist in the graph, and so could not be removed.
"""
function remove_link!(sg::SequenceDistanceGraph, source::NodeID, dest::NodeID)
    slinks = links(sg, source)
    slinkslen = length(slinks)
    filter!(!isequal(SequenceGraphLink(source, dest, 0)), slinks)
    dlinks = links(sg, dest)
    dlinkslen = length(dlinks)
    filter!(!isequal(SequenceGraphLink(dest, source, 0)), dlinks)
    return slinkslen != length(slinks) || dlinkslen != length(dlinks)
end

remove_link!(sg::SequenceDistanceGraph, lnk::DistanceGraphLink) = remove_link!(sg, source(lnk), destination(lnk))

"""
Removes all the links in the collection from and to a given nodeID.
"""
function disconnect_node!(sg::SequenceDistanceGraph, n::NodeID)
    for flink in forward_links(sg, n)
        remove_link!(sg, flink)
    end
    for rlink in backward_links(sg, n)
        remove_link!(sg, rlink)
    end
end

###
### Graph traversal
###

"""
    forward_links(sg::SequenceDistanceGraph, n::NodeID)

Get a vector of the links leaving `n` forward from that node.
"""
function forward_links(sg::SequenceDistanceGraph, n::NodeID)
    r = Vector{DistanceGraphLink}()
    nodelinks = links(sg, n)
    sizehint!(r, length(nodelinks))
    for link in nodelinks
        if is_forward_from(link, n)
            push!(r, link)
        end
    end
    return r
end

"""
    backward_links(sg::SequenceDistanceGraph, n::NodeID)

Get a vector of the links leaving `n` backwards from that node.
"""
backward_links(sg::SequenceDistanceGraph, n::NodeID) = forward_links(sg, -n)


"""
    get_next_nodes(sg::SequenceDistanceGraph, n::NodeID)

Find node IDs for forward nodes of `n`.
"""
function get_next_nodes(sg::SequenceDistanceGraph, n::NodeID)
    r = Vector{NodeID}()
    nodelinks = links(sg, n)
    sizehint!(r, length(nodelinks))
    for link in nodelinks
        if is_forwards_from(link, n)
            push!(r, destination(link))
        end
    end
    return r
end

"""
    get_previous_nodes(sg::SequenceDistanceGraph, n::NodeID)

Find node IDs for backward nodes of `n`.
"""
function get_previous_nodes(sg::SequenceDistanceGraph, n::NodeID)
    r = Vector{NodeID}()
    nodelinks = links(sg, n)
    sizehint!(r, length(nodelinks))
    for link in nodelinks
        if is_backwards_from(link, n)
            push!(r, -destination(link))
        end
    end
end

function dump_to_gfa1(sg, filename)
    fasta_filename = "$filename.fasta"
    gfa = open("$filename.gfa", "w")
    fasta = open(FASTA.Writer, fasta_filename)
    
    println(gfa, "H\tVN:Z:1.0")
    
    for nid in eachindex(nodes(sg))
        n = node(sg, nid)
        if n.deleted
            continue
        end
        println(gfa, "S\tseq", nid, "\tLN:i:", length(n.seq), "\tUR:Z:", fasta_filename)
    end
    
    for ls in links(sg)
        for l in ls
            if source(l) <= destination(l)
                print(gfa, "L\t")
                if source(l) > 0
                    print(gfa, "seq", source(l), "\t-\t")
                else
                    print(gfa, "seq", -source(l), "\t+\t")
                end
                if destination(l) > 0
                    print(gfa, "seq", destination(l), "\t+\t")
                else
                    print(gfa, "seq", -destination(l), "\t-\t")
                end
                println(gfa, distance(l) < 0 ? -distance(l) : 0, "M")
            end
        end
    end
    close(gfa)
    close(fasta)
end


    