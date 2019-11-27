## de Bruijn Graph type

A fundamental approach for de-novo gene assembly is to use make use of de Bruijn graphs.
The graph is used to represent fragments of reads (mostly starting with kmers) as vertices and
overlaps between these fragments as edges.
`DeBruijnGraph` type is a special type of `SequenceGraph`. It is also made up of two fields:

- nodes
- links

```
struct DeBruijnGraph
    nodes::Vector{SequenceGraphNode}
    links::Vector{Vector{SequenceGraphLink}}
end
```

We initialize a DeBruijnGraph using the `deBruijn_constructor` function.
This is mainly due to the fact that arbitrary links between two vertices are not allowed in the
de Bruijn graph formalism. The constructor receives as input a list of kmers and generates the deBruijn_Graph
where each kmer is a unique vertex and each overlap  of length $k-1$ is represented with an edge.

```
kmer_vector = generate_random_kmers(DNA,4,10)
dbg = deBruijn_constructor(kmer_vector)

DeBruijnGraph(SequenceGraphNode[SequenceGraphNode{Kmer{DNA,4}}(CGCC, true), SequenceGraphNode{Kmer{DNA,4}}(TCTG, true), SequenceGraphNode{Kmer{DNA,4}}(TGTG, true), SequenceGraphNode{Kmer{DNA,4}}(GAAG, true), SequenceGraphNode{Kmer{DNA,4}}(GGCA, true), SequenceGraphNode{Kmer{DNA,4}}(ACGA, true), SequenceGraphNode{Kmer{DNA,4}}(CGCT, true), SequenceGraphNode{Kmer{DNA,4}}(TCTC, true), SequenceGraphNode{Kmer{DNA,4}}(TACG, true), SequenceGraphNode{Kmer{DNA,4}}(GCAT, true)], Array{SequenceGraphLink,1}[[], [], [], [], [SequenceGraphLink(-5, 10, 1)], [], [], [], [SequenceGraphLink(-9, 6, 1)], []])
```


We can use the fastq readers to work on real data. The example below generates the set of unique kmers from 5-long reads.
The kmer size is set to 15.


```
r = FASTQ.Reader(open("URnano_ecoli.fastq", "r"))
Ecoli_reads = Vector{BioSequence{DNAAlphabet{4}}}()

## get first 5 reads
for i in 1:5
    next_seq = iterate(r)
    seq = sequence(next_seq[1])
    next_seq = iterate(r,next_seq[2])
    push!(Ecoli_reads,seq)
end
kmers = new_extract_canonical_kmers(Ecoli_reads,15)
```

The output is a set of kmers in their canonical form:

Set(Kmer{DNA,15}[AAAATATCTCGTTTT, TAAACCAGTCGCCGC, TTCGACATTACCCAG, CGCCCTGCCAGCAGT, TAATATTGTTCCATT, TGGTAATGGTCACAG, AAAAATTAAGCAGGA, ATATAAGTTATATCA, GCCCGATCTGTCTCC, GATTTCTCCGGGCCA  …  TGAGCGATTGCCTGA, CGGAGCAGCAGTGTC, AAAATCGTACATACC, GTCGCCTGATGCCTG, GTCAGCGAACCTTCC, CGCCGCTCACCGCCG, TGGATGAACGTTCAT, AATATGTCACAATTT, ATGCGATAGCAGGGG, GTAGAAAGCTCGTGG, CGATTGGTTTAAGAC])


Example data for checking node merging :

```
kmer_vector2 = Vector([DNAKmer{4}("ATTC"),DNAKmer{4}("TTCG"),DNAKmer{4}("TCGT"),
        DNAKmer{4}("AATC"),DNAKmer{4}("AATG"),DNAKmer{4}("CGTA"),DNAKmer{4}("CGTC")])
```


DeBruijnGraph for the above kmers:

```
DeBruijnGraph(Dict{Int64,SequenceGraphNode}(7=>SequenceGraphNode{Kmer{DNA,4}}(ACGA, true),4=>SequenceGraphNode{Kmer{DNA,4}}(CGTC, true),2=>SequenceGraphNode{Kmer{DNA,4}}(AATG, true),3=>SequenceGraphNode{Kmer{DNA,4}}(ATTC, true),5=>SequenceGraphNode{Kmer{DNA,4}}(CGTA, true),6=>SequenceGraphNode{Kmer{DNA,4}}(CGAA, true),1=>SequenceGraphNode{Kmer{DNA,4}}(AATC, true)), Dict(7=>[SequenceGraphLink(-7, 6, -3)],4=>[SequenceGraphLink(4, 7, -3)],2=>[],3=>[SequenceGraphLink(3, 1, -3), SequenceGraphLink(3, 2, -3)],5=>[SequenceGraphLink(5, 7, -3)],6=>[SequenceGraphLink(-6, -3, -3)],1=>[]), 4)
```


Then we apply node merging:

```
DeBruijnGraph(Dict{Int64,SequenceGraphNode}(7=>SequenceGraphNode{BioSequence{DNAAlphabet{4}}}(ACGAAT, true),4=>SequenceGraphNode{Kmer{DNA,4}}(CGTC, true),2=>SequenceGraphNode{Kmer{DNA,4}}(AATG, true),5=>SequenceGraphNode{Kmer{DNA,4}}(CGTA, true),1=>SequenceGraphNode{Kmer{DNA,4}}(AATC, true)), Dict(7=>[SequenceGraphLink(-7, 6, -3), SequenceGraphLink(-7, 1, -3), SequenceGraphLink(-7, 2, -3)],4=>[SequenceGraphLink(4, 7, -3)],2=>[],5=>[SequenceGraphLink(5, 7, -3)],1=>[]), 4)
```


So the nodes with nodeID 6 and 3 are collapsed into node 7 and we can see that both outgoing edges of 3 are given to node 7.



#### Saving to gfa

Below is an example code for saving a SequenceDistanceGraph to GFA file. This allows us to load and visualize the graph using graph visualization tools.

```
# A helpful function to let me run Bandage to visualize graphs.
const BANDAGE_BIN = "/Applications/Bandage.app/Contents/MacOS/Bandage"
function draw_graph(gr)
    filename = tempname()
    BioSequenceGraphs.dump_to_gfa1(gr, filename)
    run(`$BANDAGE_BIN image $filename.gfa $filename.png --height 500`)
    run(`open $filename.png`)
end
function show_graph(gr)
    filename = tempname()
    BioSequenceGraphs.dump_to_gfa1(gr, filename)
    run(`$BANDAGE_BIN load $filename.gfa`)
end
```
