

const NodeID = Int64

"""
A SequenceGraphNode represents nodes or vertices in a sequence graph.

Conceptually, a node in a sequence graph is a node which posesses a
biological sequence, such nodes may be connected to other nodes through
'Links'.

SequenceGraphNodes are directional. They have a start (or 'sink' makred '-'),
and they have an end (or 'source', marked '+').

Furthuremore, nodes are canonical or palindromic sequences, else they are reverted.
"""
struct SequenceGraphNode{S <: BioSequence}
    sequence::S
    active::Bool
end

sequence(sn::SequenceGraphNode) = sn.sequence

isactive(sn::SequenceGraphNode) = sn.active

function reverse_complement!(sn::SequenceGraphNode)
    reverse_complement!(sequence(sn))
end

function reverse_complement(sn::SequenceGraphNode{S}) where S <: BioSequence
    return SeqeunceNode{S}(reverse_complement(sequence(sn)), isactive(sn))
end





