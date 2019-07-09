
module GFA1

import Automa.RegExp
const re = RegExp

include("machines.jl")

end # module

function makevis(x, name)
    m = Automa.compile(x)
    write("$name.dot", Automa.machine2dot(m))
    run(`dot -Tpng -o $name.png $name.dot`)
end




cigar = re"([0-9]+[MIDNS])+"
overlap = re.cat(num, re.opt(re.cat(':', num))) | re.cat(num, ':') | re.cat(':', num) | cigar



link = re.cat('L', hspace, id, hspace, sgn, hspace, id, hspace, sgn, hspace, overlap, re.rep(re.cat(hspace, tag)))



segment_actions = Dict(
    :mark => :(@mark),
    :pos => :(pos = @relpos(p)),
    :countline => :(linenum += 1),
    :id => :(record.id = pos:@relpos(p) - 1)
    :segment => quote
        found = true
        @escape
    end
)

mutable struct SegmentRecord
    # data and filled range
    data::Vector{UInt8}
    filled::UnitRange{UInt}
    # indexes
    identifier::UnitRange{UInt}
    sequence::UnitRange{UInt}
    tags::Vector{UnitRange{UInt}}
end

function initialize!(sr::SegmentRecord)
    empty!(sr.data)
    filled = 1:0
    identifier = 1:0
    description = 1:0
    empty!(tags)
end

initcode = quote
    pos = 0
    found = false
    initialize!(record)
    cs = linenum = state
end


context = Automa.CodeGenContext()
@eval function parse_segment(g::SequenceGraph, data)
    $(Automa.generate_init_code(context, segment))
    p_end = p_eof = lastindex(data)
    $(Automa.generate_exec_code(context, machine, actions))

end
