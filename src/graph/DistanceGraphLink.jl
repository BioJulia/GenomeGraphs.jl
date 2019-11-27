###
### Link/Edge type for SequenceDistanceGraph
###

"""
Represents a single distance between two sequences in a SequenceDistanceGraph.
"""
struct DistanceGraphLink
    source::NodeID
    destination::NodeID
    dist::Int64
end

DistanceGraphLink(src::NodeID, dst::NodeID) = DistanceGraphLink(src, dst, 0)

source(l::DistanceGraphLink) = l.source
destination(l::DistanceGraphLink) = l.destination
distance(l::DistanceGraphLink) = l.dist

function Base.:(==)(x::DistanceGraphLink, y::DistanceGraphLink)
    return source(x) == source(y) && destination(x) == destination(y)
end

"Test if link `l` is a forward link leaving node `n`."
is_forwards_from(l::DistanceGraphLink, n::NodeID) = source(l) == -n
is_backwards_from(l::DistanceGraphLink, n::NodeID) = source(l) == n
