abstract type GraphEditorOp end

function index end
function input_nodes end
function input_ends end
function consumed_nodes end
function consumed_ends end

struct LinkDeletion <: GraphEditorOp
    index::Int
    src::NodeID
    dst::NodeID
end
input_ends(ld::LinkDeletion) = (src, dst)
consumed_ends(ld::LinkDeletion) = input_ends(ld::LinkDeletion)
input_nodes(ld::LinkDeletion) = ()
consumed_nodes(ld::LinkDeletion) = ()

struct NodeDeletion <: GraphEditorOp
    index::Int
    input::NodeID
    connected_ends::Vector{NodeID}
end
input_nodes(nd::NodeDeletion) = input
consumed_nodes(nd::NodeDeletion) = input_nodes(nd)
consumed_ends(nd::NodeDeletion) = nd.connected_ends

struct NodeExpansion <: GraphEditorOp
    index::Int
    node::NodeID
    input_ends::Vector{NodeID}
end
input_nodes(ne::NodeExpansion) = ne.node
consumed_nodes(ne::NodeExpansion) = input_nodes(ne)
input_ends(ne::NodeExpansion) = ne.input_ends
consumed_ends(ne::NodeExpansion) = input_ends(ne)

struct GraphEditor{G<:SequenceDistanceGraph}
    graph::G
    queued_nodes::BitVec
    queued_plus_ends::BitVec
    queued_minus_ends::BitVec
    link_deletion_queue::Vector{LinkDeletion}
    node_deletion_queue::Vector{NodeDeletion}
    next_op = Ref{Int}(0)
end

function queue_allows(ge::GraphEditor, op::GraphEditorOp)
    qn = ge.queued_nodes
    qp = ge.queued_plus_ends
    qm = ge.queued_minus_ends
    for node in input_nodes(op)
        if qn[abs(node)] || qp[abs(node)] || qm[abs(node)]
            return false
        end
    end
    for n in input_ends(op)
        qn[abs(n)] && return false
        (n < 0 && qm[-n]) && return false
        (n > 0 && qp[n]) && return false
    end
    return true
end

function queue_mark_inputs!(ge, op::GraphEditorOp)
    qn = ge.queued_nodes
    qp = ge.queued_plus_ends
    qm = ge.queued_minus_ends
    for n in consumed_nodes(op)
        qn[abs(n)] = true
        qp[abs(n)] = true
        qm[abs(n)] = true
    end
    for e in consumed_ends(op)
        if e < 0
            qm[-e] = true
        end
        if e > 0
            qp[e] = true
        end
    end
end

function queue_link_deletion!(ge::GraphEditor, src::NodeID, dest::NodeID)
    op = LinkDeletion(ge.next_op[], src, dest)
    if !queue_allows(ge, op)
        return false
    end
    ge.next_op[] += 1
    queue_mark_inputs!(ge, op)
    push!(ge.link_deletion_queue, op)
    return true
end

n_links(graph, node) = length(Graph.linksof(graph, node))

"""
    queue_node_deletion!(ge::GraphEditor, node::NodeID)

Request a GraphEditor schedules a graph node for deletion.

Returns `true` if the node deletion is compatible with previously queued operations
and has been successfully queued. Otherwise, returns `false`.
"""
function queue_node_deletion!(ge::GraphEditor, node::NodeID)
    links = Graph.linksof(ge.graph, node)
    connected_ends = [destination(l) for l in links]
    op = NodeDeletion(ge.next_op[], node, connected_ends)
    if !queue_allows(ge, op)
        return false
    end
    ge.next_op[] += 1
    queue_mark_inputs!(ge, op)
    push!(ge.node_deletion_queue, op)
    return true
end

