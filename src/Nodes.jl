

const NodeID = Int64

"""
A SequenceGraphNode represents nodes or vertices in a sequence graph.

Conceptually, a node in a sequence graph is a node which posesses a
biological sequence, such nodes may be connected to other nodes through
'Links'.

SequenceGraphNodes are directional. They have a start (marked '+'),
and they have an end (marked '-').

Furthuremore, nodes are canonical or palindromic sequences, else they are reverted.
"""
struct SequenceGraphNode{S <: BioSequence}
    sequence::S
    active::Bool
end

sequence(sn::SequenceGraphNode) = sn.sequence

isactive(sn::SequenceGraphNode) = sn.active

reverse_complement!(sn::SequenceGraphNode) = reverse_complement!(sequence(sn))

length(sn::SequenceGraphNode) = length(sequence(sn))

function iscanonical(seq::BioSequence{DNAAlphabet{2}})
    i = 1
    j = length(seq)
    @inbounds while i < j
        f = seq[i]
        r = complement(seq[j])
        f < r && return true
        r < f && return false
        i += 1
        j -= 1
    end
    return true
end

iscanonical(sn::SequenceGraphNode) = iscanonical(sequence(sn))