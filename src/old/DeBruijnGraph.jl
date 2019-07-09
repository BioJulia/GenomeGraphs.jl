"""
    DeBruijnGraph is also SequenceGraph with a special design concept
    k denotes the initial kmer length that is used during constructing the dbg
"""

struct DeBruijnGraph
    nodes::Dict{Int64,SequenceGraphNode}
    links::Dict{Int64,Vector{SequenceGraphLink}}
    k ::Int
end


nodes(dbg::DeBruijnGraph) = dbg.nodes
node(dbg::DeBruijnGraph, i::NodeID) = nodes(dbg)[abs(i)]
k_value(dbg::DeBruijnGraph) = dbg.k
links(dbg::DeBruijnGraph) = dbg.links
indegree(dbg::DeBruijnGraph,i::NodeID) = count_indegree(dbg,i)
outdegree(dbg::DeBruijnGraph,i::NodeID) = count_outdegree(dbg,i)
bothways_indegree(dbg::DeBruijnGraph,i::NodeID) = count_indegree(dbg,i)+ count_indegree(dbg,-i)
bothways_outdegree(dbg::DeBruijnGraph,i::NodeID) = count_outdegree(dbg,i)+ count_outdegree(dbg,-i)


function canonical_bio(x::BioSequence{A})where{A}
    y = reverse_complement(x)
    return x < y ? x : y
end

"""
    links(sg::SequenceGraph, node::NodeID)

Get all of the links of a Node of a sequence graph.
"""
function links(dbg::DeBruijnGraph, node::NodeID)
    l = links(dbg)
    @inbounds return l[abs(node)]
end

"""

For now this is empty. We have to decide about the design concepts first!!!
"""
function add_link!(dbg::DeBruijnGraph,l::SequenceGraphLink)

end


"""
    is_overlap(sn1::SequenceGraphNode,sn2::SequenceGraphNdoe,k::Int)

k is given as an additional information, as the length of the node labels are subject to changes during node merging
Returns true if suffix to prefix overlap of length k-1 exists
Checks overlap of k-1 long suffix of sn1 with prefix of sn2
"""
function is_overlap(sn1::SequenceGraphNode,sn2::SequenceGraphNode,k::Int)
    l1 = length(sn1)
    l2 = length(sn2)
    String(sn1.sequence)[l1-k+2:end]==String(sn2.sequence)[1:k-1]
end

"""
    find_overlaps(X::Vector{SequenceGraphNode})

Finds the overlaps of length length(Node.sequence)-1 between all Nodes in the vector X

Returns a vector of Tuples where each Tuple (i,j) corresonds to a directed edge (forward link) from i to j

For constructing the De Bruijn Graph links we need the overlapping sequence information
We first formulate the forward links as tuples with indices
Then using these tuple indices we will create the SequenceGraphLinks and finalize the Graph Construction
"""
function find_overlaps(X::Vector{SequenceGraphNode},k::Int)
    overlaps = Vector{Tuple{Int64,Int64}}()
    for i in eachindex(X)
        for j in eachindex(X)
            if is_overlap(X[i],X[j],k)
            #if String(X[i].sequence)[2:end]==String(X[j].sequence)[1:end-1]
                #println("AAAAAAAAAAAAA")
                #println("Overlap between $(X[i].sequence) and $(X[j].sequence)")
                push!(overlaps,(i,j))
            end
        end
    end
    return overlaps
end

"""
    find_overlaps(dbg::DeBruijnGraph,n::SequenceGraphNode)

Returns a list of tuples which will be used to add links to the DeBruijnGraph
finds the overlaps  of length k-1 between a node and all nodes in a dbg
"""
function find_overlaps(dbg::DeBruijnGraph,n::SequenceGraphNode,k::Int)
    nodeid = Base.length(dbg.nodes)+1
    overlaps = Vector{Tuple{Int64,Int64}}()
    for i in 1:nodeid-1
        if is_overlap(dbg.nodes[i],n,k)
        #if String(dbg.nodes[i].sequence)[2:end]==String(n.sequence)[1:end-1]
            push!(overlaps,(i,nodeid))
            println("Overlap between $(dbg.nodes[i].sequence) and $(n.sequence)")
        elseif is_overlap(n,dbg.nodes[i],k)
        #elseif String(dbg.nodes[i].sequence)[1:end-1]==String(n.sequence)[2:end]
            push!(overlaps,(nodeid,i))
            println("Overlap between $(dbg.nodes[i].sequence) and $(n.sequence)")
        end
    end
    return overlaps
end


"""
    add_node(dbg::DeBruijnGraph,n::SequenceGraphNode)

Before adding the node to the graph we find the overlaps between the newly added node and then include it
Maybe we can consider checking if it already exists
For now my plan is to add each corresponding link automatically when we add a node
"""
function add_node!(dbg::DeBruijnGraph,n::SequenceGraphNode)
    overlaps  = find_overlaps(dbg,n)
    push!(dbg.nodes,n)
    len = Base.length(dbg.nodes)
    push!(dbg.links,Vector{SequenceGraphLink}()) ## push an empty vector to the links
    for overlap in overlaps
        link = SequenceGraphLink(-overlap[1],overlap[2])
        push!(dbg.links[overlap[1]],link)
        println(link)
    end
    nodeid = Base.length(dbg.nodes)
    nodeid
end



"""
    new_deBruijn_Constructor(kmer_set::Set{Kmer{T,K}})where{T,K}

Input :  a set of kmers in their canonical form

(canonical_kmer_set function under kmer.jl can be used for creating the set out of an array of kmers)

Returns a dbg consisting of a vector of SequenceGraphNodes, a vector of SequenceGraphLinks and an integer k for the kmer length

New Constructor for the DeBruijnGraph where we simultaneously consider a kmer and its reverse complement
We define the beginning of the canonical form of the kmer as the (+) end
And end of the canonical form as the (-) ends

"""
function new_deBruijn_Constructor(kmer_set::Set{Kmer{T,K}})where{T,K}
    fw_nodes = Vector{Tuple{Kmer{T,K-1},Int64}}()
    bw_nodes = Vector{Tuple{Kmer{T,K-1},Int64}}()

    nodes = Vector{SequenceGraphNode}()
    links  = Vector{Vector{SequenceGraphLink}}()

    i = 1
    flag = 0
    ## adding unique kmers to graph in their canonical form
    for kmer in kmer_set
        can_kmer = canonical(kmer)
        push!(nodes,SequenceGraphNode(can_kmer,true))
        pref = Kmer{T,K-1}(String(can_kmer)[1:end-1])
        suf = Kmer{T,K-1}(String(can_kmer)[2:end])


        ### adding both ways if pref == reverse_complement(pref)
        if canonical(pref) == pref ## add prefix to forward nodes
            push!(fw_nodes,(pref,i))
            if pref == reverse_complement(pref)
                push!(bw_nodes,(pref,i))
            end
        else
            push!(bw_nodes,(canonical(pref),i))
        end

        a = 1
        ### adding both ways if suf == reverse_complement(suf)
        if canonical(suf) == suf ## add suffix to backward nodes (outgoing edge)
            push!(bw_nodes,(suf,-i))
            if suf == reverse_complement(suf)
                push!(fw_nodes,(suf,-i))
            end
        else
            push!(fw_nodes,(canonical(suf),-i))
        end
        i= i + 1
    end
    prev = abs(last(bw_nodes[1]))
    links_ = Vector{SequenceGraphLink}()
    for kbn in bw_nodes
        if abs(last(kbn))!=prev
            push!(links,links_)
            for i in prev+1:abs(last(kbn))-1
                push!(links,Vector{SequenceGraphLink}())
            end
            links_ = Vector{SequenceGraphLink}()
        end
        for kfn in fw_nodes
            if first(kfn)==first(kbn)
                push!(links_,SequenceGraphLink(last(kbn),last(kfn),-K+1)) ## i did not quiet get -k+1
            end
        end
        prev = abs(last(kbn))
    end
    push!(links,links_)
    nodes_new = Dict{Int64,SequenceGraphNode}()
    links_new = Dict{Int64,Vector{SequenceGraphLink}}()
    for (node_id,node) in enumerate(nodes)
        nodes_new[node_id] = node
        links_new[node_id] = links[node_id]
    end
    dbg = DeBruijnGraph(nodes_new,links_new,K)
end



"""
    deBruijn_constructor(kmer_vector::Vector{Kmer{T,K}}) where{T<:NucleicAcid,K}

Returns a DeBruijnGraph constructed by the kmers

For now lets assume that we have Kmers prepared already
So for an unknown DNA sequence we have Spectrum(s,k) where Spectrum(s,l)
is the multiset of n-l+1 l-mers in s.

Our De Bruijn Graph constructor takes as input some kmers and puts a directed link between
every two nodes n1,n2 such that n1[2:] = n2[:k-1], i.e overlap(n1,n2) == k-1
first we find the overlaps between each node by doing an computationally expensive exhaustive search
over all pairs

deBruijn_Constructor takes as input some number of kmers kmer_i, generates Sequence Graph Nodes n_i
such that each n_i.sequence = kmer_i

The constructor generates  nodes in the same order in the input kmer vector

Also the constructor adds empty  vector for each node with no forward link
So the number of vectorsin Links match with number of nodes in Nodes
"""
function deBruijn_constructor(kmer_vector::Vector{Kmer{T,K}}) where{T<:NucleicAcid,K}
    Nodes = Vector{SequenceGraphNode}()
    for kmer in kmer_vector
        node = SequenceGraphNode(kmer,true)
        push!(Nodes,node)
    end
    overlaps = find_overlaps(Nodes,K)
    Links = Vector{Vector{SequenceGraphLink}}()
    prev_i = 1 ## initial NodeID
    current_node_vector = Vector{SequenceGraphLink}()
    for overlap in overlaps
        if overlap[1]!=prev_i
            push!(Links,current_node_vector)
            current_node_vector = Vector{SequenceGraphLink}()
            for i in prev_i:overlap[1]-2
                push!(Links,Vector{SequenceGraphLink}())
            end
            prev_i = overlap[1]
        end
        ## links are generated from - end of outgoing node to the + end of the incoming one
        link = SequenceGraphLink(-overlap[1],overlap[2],1) ## right now dist is initialized as 1
        push!(current_node_vector,link)
    end
    push!(Links,current_node_vector)
    for i in prev_i+1:size(Nodes)[1]
        push!(Links,Vector{SequenceGraphLink}())
    end
    deBruijn_Graph = DeBruijnGraph(Nodes,Links,K)
    deBruijn_Graph
end



# Kmer Enumeration
"""
    new_extract_canonical_kmers(read_set::Set{BioSequence{A}},k::Int64)where{A}

    Return a set of unique kmers in their canonical form which will be used to generate the dbg from a set of read_set

    Input: A set of reads (BioSequences)

    Faster way of getting all kmers using the kmer iterators in eachkmer.jl

"""
## this is the new version of extracting canonical kmers
## we use the kmer iterators instead of manually extracting each kmer
function new_extract_canonical_kmers(read_set::Set{BioSequence{A}},k::Int64)where{A}
    T = eltype(A)
    kmer_set = Set{Kmer{T,k}}()
    for seq in read_set
        kmer_iterator = each(Kmer{T,k},seq)
        kmer_state = iterate(kmer_iterator) ## initialize
        while kmer_state!=nothing
            kmer = kmer_state[1][2]
            push!(kmer_set,kmer)
            start_ind = kmer_state[1][1]
            new_state= (start_ind,1,UInt64(0))
            kmer_state = iterate(kmer_iterator,new_state)
        end

    end
    return kmer_set
end

## Error Correction !!!
# Very important step for contig formalization


"""
    get_parents(node_id::Int64,dbg::DeBruijnGraph)

returns the list of node_ids who has outgoing edges to the query nodes

Exhaustive search over all links on the DBG
We need a faster way of implementing this query
"""

function get_parents(node_id::Int64,dbg::DeBruijnGraph)
    parents = Vector{Int64}()
    for link_list in links(dbg)
        for link in link_list[2]
            if abs(destination(link))==node_id
                push!(parents,abs(link_list[1]))
            end
        end
    end
    return parents
end


"""
    delete_tips(dbg::DeBruijnGraph)

deletes the dead-end tips in a dbg

dead-end tip is defined as tips that have no outgoing edge and a single incoming edge.
The source node of the incoming edge has multiple outgoing edges (hopefully two)
candidates store all dead-end nodes.
We check whether the other outgoing edge of the source node is also dead-end. We delete if it is not.
"""
function delete_tips(dbg::DeBruijnGraph)
    candidates = Vector{Int64}()
    for link in links(dbg)
        if Base.length(link[2])==0
            push!(candidates,link[1])
        end
    end
    println("Candidates")
    println(candidates)
    to_be_removed = Vector{Int64}()
    for cand in candidates
        parents = get_parents(cand,dbg)
        for parent in parents
            if Base.length(links(dbg)[parent])>1
                for link in links(dbg)[parent]
                    if abs(destination(link))!=cand && !(abs(destination(link)) in candidates)
                        push!(to_be_removed,cand)
                    end
                end
            end
        end
    end
    println("To be removed")
    println(to_be_removed)
    for i in to_be_removed
        delete!(nodes(dbg2),abs(i))
        delete!(links(dbg2),abs(i))
    end
end


# Kmer Enumeration
"""
    extract_canonical_kmers(read_set::Set{BioSequence{A}},k::Int64)where{A}

    Return a set of unique kmers in their canonical form which will be used to generate the dbg from a set of read_set

    Input: A set of reads (BioSequences)

    Very slow way of getting all kmers

"""
function extract_canonical_kmers(read_set::Set{BioSequence{A}},k::Int64)where{A}
    T = Base.eltype(A)
    kmer_set = Set{Kmer{T,k}}()
    for seq in read_set
        start = 1
        fin = start+k-1
        while fin <= Base.length(seq)
            kmer = canonical(Kmer{T,k}(BioSequence{A}(seq,start:fin)))
            push!(kmer_set,kmer)
            start +=1
            fin +=1
        end
    end
    return kmer_set
end

# Path merging

"""
    merge_simple_paths(dbg,simple_paths;alp=DNAAlphabet{4})

Returns the updated dbg

Gets as input a dbg and a list of lists where each list contain > 0 nodeID's representing the nodes on the simple path
Simply concatenate each sequence on the path to a single sequence
Delete old intermediate links and transfer the outgoing edges of the end node if any.

"""
function merge_simple_paths(dbg::DeBruijnGraph,simple_paths;alp=DNAAlphabet{4})
    overlap = dbg.k - 1
    nodes_ = nodes(dbg)
    links_ = links(dbg)
    for path in simple_paths## list of simple paths as nodeids (with direction)
        last_ind = path[1]
        ## find the new node label
        if Base.length(path)==1
            continue
        end
        if path[1] > 0
            seq = reverse_complement(sequence(nodes_[abs(path[1])]))
        else
            seq = sequence(nodes_[abs(path[1])])
        end
        new_seq = BioSequence{alp}(seq)

        ## remove outgoing edge from the start node to next path
        next_node = path[2]
        link_index = find_link_index(links_[abs(path[1])],next_node) ## can be used if we force removal of a non-simple path
        deleteat!(links_[abs(path[1])],1)

        for nodeid in path[2:end-1]
            if nodeid < 0
                seq = reverse_complement(sequence(nodes_[abs(nodeid)]))
            else
                seq = sequence(nodes_[abs(nodeid)])
            end
            seq = sub_seq(seq,overlap+1)
            new_seq = BioSequence{alp}(new_seq, seq)
            delete!(nodes_,abs(nodeid))
            delete!(links_,abs(nodeid))
            last_ind = nodeid
        end
        last_ind = path[end]##keep the links of the last node
        if last_ind < 0
            seq = reverse_complement(sequence(nodes_[abs(last_ind)]))
        else
            seq = sequence(nodes_[abs(last_ind)])
        end
        seq = sub_seq(seq,overlap+1)
        new_seq = BioSequence{alp}(new_seq, seq)
        nodes_[abs(path[1])] = SequenceGraphNode(new_seq,true)
        delete!(nodes_,abs(last_ind))
        ## find the new links
        for link in links(dbg)[abs(last_ind)]
            push!(links_[abs(path[1])],SequenceGraphLink(path[1],destination(link),distance(link)))
        end
        delete!(links_,abs(last_ind))
    end
    dbg
end

# Simple path finders for unitigging
# -----

"""
    simple_path_finder(dbg)

Returns a list of paths which correspond to unitigs that can not be further extended

More testing is necessary

The main idea is that there  can be  at most 1  simple path originating from each node with outdegree 1

So an exhaustive search is not very expensive if we only look at nodes with incoming > 1 or  incoming == 0 and out==1
Looking at nodes with incoming >  1 will give us the maximal unitigs
On the other hand, looking at nodes with outgoing = 1 will give all simple paths
"""

function simple_path_finder(dbg::DeBruijnGraph)
    println("ADRARAR")
    simple_path = Vector()
    links1 = links(dbg)
    for node_id in 1:Base.length(nodes(dbg))
        path = Vector()
        parent = node_id
        flag = true

        ## start of a maximal unitig in_degree >2 or 0 and out_degree ==1
        if (bothways_indegree(dbg,node_id)==0 || bothways_indegree(dbg,node_id)>1 )&& bothways_outdegree(dbg,node_id)==1
            sign1 = sign(source(links1[node_id][1]))
            push!(path,sign1*parent)
            parent = destination(links1[node_id][1])
            while bothways_outdegree(dbg,parent)==1 && bothways_indegree(dbg,parent)==1
                sign1 = sign(source(links1[abs(parent)][1]))
                push!(path,parent)
                parent = destination(links1[abs(parent)][1])
            end

            ## if the indegree==1 extend the path to include the final kmer at then end of the path
            ## o/w skip
            if bothways_indegree(dbg,parent)==1
                push!(path,parent)
            end
            println("Maximal path found!")
            println(path)
            push!(simple_path,path)
        end
    end
    simple_path
end


# Query Functions
# --------

# Counters
# ------
# Counters count the indegree and outdegree separately for the (+) and  (-) nodes
# edges to and out of the  (+) end  represent the path through canonical kmer
# edges to and out of the  (-) end represent the path through non-canonical kmer

"""
    count_indegree(dbg::DeBruijnGraph,n::NodeID)

Returns the number of incoming edges to the SequenceGraphNode with id NodeID on the DeBruijnGraph

Makes an exhaustive search on all the SequenceGraphLinks on the dbg
Now checking only the source end of the node we have to discuss about the design
Maybe we can store incoming and outgoing edges for each node for O(N) query time where N denotes number of SequenceGraphNodes
"""
function count_indegree(dbg::DeBruijnGraph,n::NodeID)
    dest = n
    in_degree = 0
    for i in 1:Base.length(nodes(dbg))
        links_ = links(dbg,i)
        for l in links_
            if destination(l)==dest  ## checking only the source end of the node, not sure!
                in_degree +=1
            end
        end
    end
    in_degree
end



"""
    count_outdegree(dbg::DeBruijnGraph,n::NodeID)

Returns the number of incoming edges to the SequenceGraphNode with id NodeID on the DeBruijnGraph

Makes an exhaustive search on all the SequenceGraphLinks on the dbg
Now checking only the sink end of the node we have to discuss about the design
O(K) query time
"""
function count_outdegree(dbg::DeBruijnGraph,n::NodeID)
    source_ =  n  ## checking the sink end of the node for outgoing edges
    out_degree = 0
    for link in links(dbg,n)
        if source_ == abs(source(link))
            out_degree +=1
        end
    end
    out_degree
end


# Path queries
# ------


"""

Returns if a given sequence is in a dbg  and also all the internal nodes have in_degree and  out_degree equal to one.


"""

function is_simple_path(seq::Sequence,dbg::DeBruijnGraph)
    is_path,path = is_a_path2(seq,dbg,min_match=k_value(dbg)-1)
    print(path)
    simple_path = []
    if is_path
        if bothways_outdegree(dbg,abs(path[1]))>1
            println("Out degree of node $(path[1]) is $(bothways_outdegree(dbg,abs(path[1])))")
            return false
        end
        push!(simple_path,abs(path[1]))
        for nodeid in path[2:end-1]
            indeg = bothways_indegree(dbg,abs(nodeid))
            outdeg = bothways_outdegree(dbg,abs(nodeid))
            if indeg==1 &&outdeg==1
                push!(simple_path,abs(nodeid))
            else
                println("Node $(nodeid) has $(indeg) indegrees and $(outdeg) outdegreess")
                return false,simple_path

            end

        end
        if Base.length(path)==1
            return true,simple_path
        end
        if bothways_indegree(dbg,abs(path[end]))==1
            push!(simple_path,abs(path[end]))
            return true,simple_path
        else
            return false,simple_path
        end
    end
    return false
end





"""
    is_a_path2(seq::Sequence,dbg::DeBruijnGraph;min_match=3)

Version two for is_a_path for the new dbg designed
Soon will replace version 1

Returns the NodeIDs of all matches and returns true if the path exists completely in dbg

If complete match is not found still returns the indexes of all the matches

min_match is the initial kmer length k used during dbg construction
After node merging we still have k-1 overlaps between nodes because we only merge the simple paths.
Start from a node and traverse its links.
Stop at the first false during traversal as there can not be another path starting from another node!


Making use of is_suffix from Nodes.jl and is_match from sequence.jl
TO-DOs :
    1 - Make sure these dependencies function properly
    2 - Check if a kmer  is contained completely in a node
"""
function is_a_path2(seq::Sequence,dbg::DeBruijnGraph;min_match=3)
    #seq = canonical_bio(seq)
    println(seq)
    index = 0
    indexes = Vector{Int64}()
    nodes_ = nodes(dbg)
    match = -1
    for i in 1:Base.length(nodes_)
        match =  is_suffix2(seq,nodes_[i],direction=1,min_match=min_match)
        match_rev = is_suffix2(seq,nodes_[i],direction=-1,min_match=min_match)
        if match == Base.length(seq) || match_rev == Base.length(seq)
            push!(indexes,i)
            return true, indexes
        end

        if match!=-1
            index = -i
            println(nodes_[i])
            break
        elseif match_rev!=-1
            index = i
            match = match_rev
            seq1 = reverse_complement(sequence(nodes_[i]))
            println("Reverse complement Match: "*String(seq1))
            break
        end
    end
    if index==0
        return false

    else ## traverse on the children until remaining is smaller than or equal to min_match length
        ## since kmers overlap with k-1 min_match+1 overlap is carried to the next step
        push!(indexes,index)
        remaining_seq = sub_seq(seq,match-min_match+1)
        rem_len = Base.length(remaining_seq)
        println("Remaining: "*String(remaining_seq))
        println("Remaining length : "* string(Base.length(remaining_seq)))
        println("Index: $(index)")
        #index = index * -1
        while rem_len > 0 ## the remaining must match the next node from index 1
            found = false
            for link in links(dbg)[abs(index)]
                if source(link)==index ## if the direction  is correct
                    dest_ind  = destination(link)
                    node = nodes_[abs(dest_ind)]
                    println("Destination $(dest_ind)")
                    if dest_ind > 0
                        node_seq = sequence(node)
                    else
                        node_seq = reverse_complement(sequence(node))
                    end
                    println("Node sequence : $(node_seq)")
                    min_length = min(rem_len,Base.length(node_seq))
                    if is_match(remaining_seq,1,node_seq,1,min_length)
                        found = true
                        push!(indexes, dest_ind)
                        index = -dest_ind
                        println(nodes_[abs(index)])
                        if min_length ==rem_len
                            println("Match Found!!")
                            return true,indexes
                        end
                        remaining_seq = sub_seq(remaining_seq,min_length-min_match+1)
                        rem_len = Base.length(remaining_seq)
                        println("Remaining : $(remaining_seq)")
                        continue
                    end
                end
            end
            if found ==false
                return false,indexes
            end
        end
    end
    return true,indexes
end

"""
    is_a_path(seq::Sequence,dbg::DeBruijnGraph;min_match=3)

Returns the NodeIDs of all matches and returns true if the path exists completely in dbg

If complete match is not found still returns the indexes of all the matches

min_match is the initial kmer length k used during dbg construction
After node merging we still have k-1 overlaps between nodes because we only merge the simple paths.
Start from a node and traverse its links.
Stop at the first false during traversal as there can not be another path starting from another node!


Making use of is_suffix from Nodes.jl and is_match from sequence.jl
TO-DOs :
    1 - Make sure these dependencies function properly
    2 - Check if a kmer  is contained completely in a node
"""
function is_a_path(seq::Sequence,dbg::DeBruijnGraph;min_match=3)
    index = -1
    indexes = Vector{Int64}()
    nodes_ = nodes(dbg)
    match = -1
    for i in 1:Base.length(nodes_)
        match =  is_suffix(seq,nodes_[i],min_match =min_match)
        if match!=-1
            index = i
            println(nodes_[i])
            break
        end
    end
    if index==-1
        return false
    else ## traverse on the children until remaining is smaller than or equal to min_match length
        ## since kmers overlap with k-1 min_match+1 overlap is carried to the next step
        push!(indexes,index)
        remaining_seq = sub_seq(seq,match+1-min_match+1)
        rem_len = Base.length(remaining_seq)
        println("Remaining: "*String(remaining_seq))
        println("Remaining length : "* string(Base.length(remaining_seq)))
        while rem_len > 0 ## the remaining must match the next node from index 1
            found = false
            for link in links(dbg)[index]
                node = nodes_[destination(link)]
                node_seq = sequence(node)
                min_length = min(rem_len,Base.length(node_seq))
                println(min_length)
                if is_match(remaining_seq,1,node_seq,1,min_length)

                    found = true
                    index = destination(link)
                    push!(indexes,index)
                    println(nodes_[index])
                    if min_length ==rem_len
                        println("Match Found!!")
                        return true,indexes
                    end
                    remaining_seq = sub_seq(remaining_seq,rem_len-min_match+1)
                    rem_len = Base.length(remaining_seq)
                    println(remaining_seq)
                    continue
                end
            end
            if found ==false
                return false,indexes
            end
        end
    end
    return true,indexes
end
