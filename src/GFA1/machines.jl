
function makevis(x, name)
    m = Automa.compile(x)
    write("$name.dot", Automa.machine2dot(m))
    run(`dot -Tpng -o $name.png $name.dot`)
end

function appendfrom!(dst, dpos, src, spos, n)
    if length(dst) < dpos + n - 1
        resize!(dst, dpos + n - 1)
    end
    copyto!(dst, dpos, src, spos, n)
    return dst
end

newline = let
    lf = re" "
    lf.actions[:enter] = [:countline]

    re.cat(re.opt('\r'), lf)
end

segment = let
    hspace = re"[ \t\v]"

    seqid = re"[!-~]+"
    seqid.actions[:enter] = [:pos]
    seqid.actions[:exit] = [:seg_id]

    letters = re"[A-Za-z=.]+"
    letters.actions[:enter] = [:filledseq]

    emptyseq = re"\*"
    emptyseq.actions[:enter] = [:emptyseqs]

    sequence = letters | emptyseq
    sequence.actions[:enter] = [:pos]
    sequence.actions[:exit] = [:seg_seq]

    tag = let
        key = re"[A-Za-z][A-Za-z0-9_]"
        val = re.alt(
            re"A:[!-~]",
            re"i:[-+]?[0-9]+",
            re"f:[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?",
            re"Z:[ !-~]*",
            re"H:([0-9A-F][0-9A-F])*",
            re"B:[cCsSiIf](,[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)+")

        re.cat(key, ':', val)
    end
    tag.actions[:enter] = [:pos]
    tag.actions[:exit] = [:seg_tags]

    re.cat('S', hspace, seqid, hspace, sequence, re.rep(re.cat(hspace, tag)))
end

segment.actions[:enter] = [:mark]
segment.actions[:exit] = [:segment]

file = re.rep(re.cat(segment, newline))

segment_actions = Dict(
    :mark => :(@mark),
    :pos => :(pos = @relpos(p)),
    :countline => :(linenum += 1),
    :seg_id => :(segname = pos:@relpos(p - 1)),
    :seg_seq => :(seq = pos:@relpos(p - 1)),
    :emptyseq => :(println("Emptysequence!")),
    :seg_tags => :(tags = pos:@relpos(p - 1)),
    :segment => quote
        println("Segment!")
        appendfrom!(recorddata, 1, data, @markpos, p-@markpos)
        println(recorddata[segname])
        println(recorddata[seq])
        println(recorddata[tags])
    end
)

context = Automa.CodeGenContext(generator=:goto)
Automa.Stream.generate_reader(
    :gosegments,
    machine,
    actions = segment_actions,
    context = context
) |> eval













link = let
    hspace = re"[ \t\v]"

    fromid = re"[!-~]+"
    fromid.actions[:enter] = [:pos]
    fromid.actions[:exit] = [:from_id]

    fromsgn = re"[-+]"
    fromsgn.actions[:enter] = [:pos]
    fromsgn.exit[:exit] = [:from_sgn]

    toid = re"[!-~]+"
    toid.actions[:enter] = [:pos]
    toid.actions[:exit] = [:to_id]

    tosgn = re"[-+]"
    tosgn.actions[:enter] = [:pos]
    tosgn.actions[:exit] = [:to_sgn]


end

link.actions[:enter] = [:mark]
link.actions[:exit] = [:link]

macro re_str(line)
    return line
end
