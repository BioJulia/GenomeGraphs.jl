

"The NodeView provides a read-only comfortable interface for graph traversal."
struct NodeView{G<:SequenceDistanceGraph}
    ws::WorkSpace
    id::NodeID
    graph::G
end

"""
    Get a view of node `nodeid` in the base SequenceDistanceGraph of the
    WorkSpace `ws`.
"""
function node(ws::WorkSpace, nodeid::NodeID)
    return NodeView(ws, nodeid, ws.sdg)
end