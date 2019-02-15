
function load_from_gfa(filename::String)
    
    graph = SequenceGraph(SequenceGraphNode[], SequenceGraphLink[])
    oldnames_to_ids = Dict{String, NodeID}()
    oldnames = Vector{String}()
    
    gfaf = open(filename, "r")
    if eof(gfaf)
        throw(ArgumentError("GFA file is empty"))
    end
    line = readline(gfaf)
    if line != "H\tVN:Z:1.0"
        @warn "The first line of the GFA does not correspond to GFA1"
    end
    faf = open(FASTA.Reader, string(splitext(filename), ".fasta"))
    if eof(faf)
        throw(ArgumentError("FASTA file is empty"))
    end
    
    # Load all sequences from fasta, if they're not canonical, flip them and
    # remember they are flipped.
    
    rcnodes = 0
    
    farec = FASTA.Record()
    while !eof(faf)
        read!(faf, farec)
        nodename = FASTA.identifier(farec)
        nodeseq = FASTA.sequence(BioSequence{DNAAlphabet{2}}, farec)
        if haskey(oldnames_to_ids, nodename)
            error("sequence ", nodename, " is already defined.")
        end
        oldnames_to_ids[nodename] = add_node!(graph, Node(seq, true));
        push!(oldnames, nodename)
        # Canonicalize the sequence
        if !iscanonical(nodeseq)
            reverse_complement!(nodeseq)
            oldnames_to_ids[nodename] = -oldnames_to_ids[nodename]
            rcnodes += 1
        end
    end
    @info string(length(nodes(graph)), " nodes loaded (", rcnodes, " canonised).")
    
    dist_egt0 = 0
    lcount = 0
    # Load and store all connections
    for line in eachline(gfaf)
        tokens = split(line)
        recordtype = tokens[1]
        
        if recordtype == "S"
            id = tokens[2]
            seq = tokens[3]
            len = tokens[4]
            if seq != '*'
                error("sequences should be in a seperate fasta file.")
            end
            # Check equal length seq and node length reported in gfa
            if haskey(oldnames_to_ids, id)
                if parse(Int, len[6:end]) != length(node(graph, oldnames_to_ids[id]))
                    error(string("Different length in node and fasta for sequence: ", id, " -> gfa: ",  parse(Int, len[6:end]),  ", fasta: ", length(node(graph, oldnames_to_ids[id]))))
                end
            end
        elseif recordtype == "L"
            source = tokens[2]
            sdir = tokens[3]
            dest = tokens[4]
            ddir = tokens[5]
            cigar = tokens[6]
            if haskey(oldnames_to_ids, source)
                oldnames_to_ids[source] = add_node!(graph, Node(BioSequences{DNAAlphabet{2}}()))
                @info string("added source ", source)
            end
            if haskey(oldnames_to_ids, dest)
                oldnames_to_ids[dest] = add_node!(graph, Node(BioSequences{DNAAlphabet{2}}()))
                @info string("added dest ", dest)
            end
            source_id = oldnames_to_ids[source]
            dest_id = oldnames_to_ids[dest]
            # Go from GFA's "Links as paths" to a normal "nodes have sinks (-/start/left) and sources (+/end/right)"
            if sdir == "+"
                source_id = -source_id
            end
            if ddir == "-"
                dest_id = -dest_id
            end
            dist = 0
            if length(cigar) > 1 && cigar[end] == 'M'
                # TODO: better checks for M
                dist = -parse(Int, cigar[1:end-1])
                if dist >= 0
                    dist_egt0 += 1
                end
            end
            add_link!(graph, source_id, dest_id, dist)
            lcount += 1
        end
    end
    close(gfaf)
    if dist_egt0 > (lcount / 2)
        @warn string("Warning: The loaded graph contains ", dist_egt0, " non-overlapping links out of ", lcount)
    end
    @info string(length(nodes(graph)) - 1, " nodes after connecting with ", lcount, "links.")
end