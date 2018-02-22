
"""
A SequenceGraphLink represents a link between two nodes
(which are represented by the `SequenceGraphNode` type) in a SequenceGraph.

Conceptually, links are basically a pair of nodes with signs.
A positive node ID (e.g. 5 or +5) denotes the connection uses a node's source (end).
Conversely, a negative node ID (e.g - 5) denotes the conection uses a node's sink (start).
"""
struct SequenceGraphLink
    source::Int64
    destination::Int64
    dist::Int64
end



