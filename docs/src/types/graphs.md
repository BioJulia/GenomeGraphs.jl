# The Sequence Distance Graph (SDG)

```@docs
SequenceDistanceGraph
```

# Querying an SDG for basic properties

```@docs
nodes
n_nodes
each_node_id
links
node
sequence
```

# Manually editing an SDG by manipulating nodes and links

It is not recommended you do this if you are a high level user.
However these small editing operations are required for developers.
If you find yourself needing these methods, you will have to explicitly import
them, as they are not exported from the module.

If you find yourself wanting to edit the graph manually, it's a good idea to
ask the package authors listed in the `Project.toml` or `.github/CODEOWNERS`.

```@docs
add_node!
remove_node!
add_link!
remove_link!
disconnect_node!
```

#Code example for testing the error correction and simplifying the graph

In the example below we generate a toy graph from 4-mers which contain a bubble.
As for now, we make use of randomly generated coverage information for
```
using BioSequences
using BioSequenceGraphs

kmerlist = Vector{DNAKmer{4}}([DNAKmer{4}("ATAC"),DNAKmer{4}("TACG"),DNAKmer{4}("TACC"),DNAKmer{4}("ACGA"),DNAKmer{4}("CGAA"),DNAKmer{4}("GAAT"),DNAKmer{4}("AATC"),DNAKmer{4}("ACCA"),DNAKmer{4}("CCAA"),DNAKmer{4}("CAAT")])
sdg = SequenceDistanceGraph(kmerlist,true)
```

The new SequenceDistanceGraph constructor takes as input a list of kmerlist (not necessarily sorted or in canonical form) and a boolean for whether to do error correction or not. Then by using the coverage information generated randomly, simplifies the graph.
The resulting SDG before and after the error correction steps are as follows:

***Before:***
```
SequenceDistanceGraph{BioSequence{DNAAlphabet{2}}}(SDGNode{BioSequence{DNAAlphabet{2}}}[SDGNode{BioSequence{DNAAlphabet{2}}}(AATC, false), SDGNode{BioSequence{DNAAlphabet{2}}}(ATAC, false), SDGNode{BioSequence{DNAAlphabet{2}}}(ATTCGTA, false), SDGNode{BioSequence{DNAAlphabet{2}}}(ATTGGTA, false)], Array{DistanceGraphLink,1}[[DistanceGraphLink(1, 3, -3), DistanceGraphLink(1, 4, -3)], [DistanceGraphLink(-2, -4, -3), DistanceGraphLink(-2, -3, -3)], [DistanceGraphLink(3, 1, -3), DistanceGraphLink(-3, -2, -3)], [DistanceGraphLink(4, 1, -3), DistanceGraphLink(-4, -2, -3)]])
```

***After:***
```
SequenceDistanceGraph{BioSequence{DNAAlphabet{2}}}(SDGNode{BioSequence{DNAAlphabet{2}}}[SDGNode{BioSequence{DNAAlphabet{2}}}(GATTCGTAC, false)], Array{DistanceGraphLink,1}[[]])
```

So the algorithm first generates a sdg with collapsing all simple paths using the initial kmer list. Then at the next stage the low covered bubble branch is removed and path collapsing is applied one more time.
This results in a sdg with a whole single node, a perfect unitig.


##TODO:

Implement the sdg constructor to take as input a list of reads instead of kmers and generate all kmer + coverage information from them.
