
struct SequenceGraphLink
    source::NodeID
    destination::NodeID
    dist::Int64
end

## enable initializing links without a distance
SequenceGraphLink(s::NodeID,d::NodeID)=SequenceGraphLink(s,d,-1)

source(l::SequenceGraphLink) = l.source
destination(l::SequenceGraphLink) = l.destination
distance(l::SequenceGraphLink) = l.dist

is_forward_link(l::SequenceGraphLink, n::NodeID) = source(l) == n

@inline function Base.isequal(a::SequenceGraphLink, b::SequenceGraphLink)
    return source(a) == source(b) && destination(a) == destination(b)
end


function find_link_index(links_vec::Vector{SequenceGraphLink},nodeid)
    for (i,link) in enumerate(links_vec)
        if destination(link) == nodeid
            return i
        end
    end
end
