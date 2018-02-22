

struct SequenceNode{S <: BioSequence}
    sequence::S
    active::Bool
end

sequence(sn::SequenceNode) = sn.sequence

isactive(sn::SequenceNode) = sn.active

function reverse_complement!(sn::SequenceNode)
    reverse_complement!(sequence(sn))
end

function reverse_complement(sn::SequenceNode{S}) where S <: BioSequence
    return SeqeunceNode{S}(reverse_complement(sequence(sn)), isactive(sn))
end





