###
### Node type for SequenceDistanceGraph
###

"""
The SDGNode type represents a node in a SequenceDistanceGraph.

At present it contains only two fields, first it holds an instance of a BioSequences
Sequence type.

Secondly, it tracks a flag which indicates if the node has been
deleted or not.

!!! note
    The deleted flag allows us to mark nodes in the graph as deleted, which can be of
    help in some algorithms where the graph structure is being edited (merging nodes for example).

    Actually deleting the node would shift node IDs and require redoing all links in the graph
    and so on.

    So just marking a node as deleted and not using it anymore is a lazy but sometimes
    helpful choice.
"""
struct SDGNode{S}
    seq::S
    deleted::Bool
end

function empty_node(::Type{S}) where {S <: Sequence}
    return SDGNode{S}(empty_seq(S), true)
end

# TODO: This is a hacked copy of the dna string literal macro from BioSequences,
# except it creates 2-bit based DNA sequences rather than 4 bit based ones.
# This ability to choose the bit encoding should make its way to BioSequences.jl
# in the future, but for now, it's here.
#
# I basically want this as it lets me create a single literal empty sequence, shared
# by all deleted SDG nodes. Rather than having each deleted SDG node create a new empty
# sequence.
macro dna2_str(seq, flag)
    if flag == "s"
        return BioSequence{DNAAlphabet{2}}(BioSequences.remove_newlines(seq))
    elseif flag == "d"
        return quote
            BioSequence{DNAAlphabet{2}}($(BioSequences.remove_newlines(seq)))
        end
    end
    error("Invalid DNA flag: '$(flag)'")
end

@inline empty_seq(::Type{BioSequence{DNAAlphabet{2}}}) = dna2""s

@inline is_deleted(n::SDGNode{S}) where {S<:Sequence} = n.deleted
@inline length(n::SDGNode{S}) where {S<:Sequence} = length(n.seq)