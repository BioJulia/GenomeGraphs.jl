

struct SequenceGraphLink
    source::NodeID
    destination::NodeID
    dist::Int64
end

source(l::SequenceGraphLink) = l.source
destination(l::SequenceGraphLink) = l.destination
distance(l::SequenceGraphLink) = l.dist

is_forward_link(l::SequenceGraphLink, n::NodeID) = source(l) == -n

@inline function Base.isequal(a::SequenceGraphLink, b::SequenceGraphLink)
    return source(a) == source(b) && destination(a) == destination(b) 
end