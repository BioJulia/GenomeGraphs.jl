# Package Guide

## Installation

GenomeGraphs is made available to install through BioJulia's package registry.

Julia by default only watches the "General" package registry, so before you
start, you should add the BioJulia package registry.

Start a julia terminal, hit the `]` key to enter pkg mode (you should see the
prompt change from `julia>` to `pkg>`), then enter the following command:

```
pkg> registry add https://github.com/BioJulia/BioJuliaRegistry.git
```

After you've added the registry, you can install GenomeGraphs from the julia REPL.
Press `]` to enter pkg mode again, and enter the following:

```
pkg> add GenomeGraphs
```

If you are interested in the cutting edge of the development, please check out
the master branch to try new features before release.

## Creating your first WorkSpace

You can create an empty genome graph workspace with the empty constructor:

```julia
using GenomeGraphs

ws = WorkSpace()
```

Obviously this workspace is quite useless on its own!

You need a genome graph and information to project onto it before you can
do any exploration or analysis. There are a few ways to get your first graph:

1. Construct a de-Bruijn graph from raw sequencing reads.
2. Load a graph (such as produced from another assembler) from a GFAv1 file

## Constructing your first de-Bruijn graph

Let's see how to do option number 1, and construct a de-Bruijn graph from
raw sequencing reads. This can be achieved with a few simple steps:

1. Prepare the sequencing reads & build a datastore.
2. Add the read datastore to a WorkSpace.
3. Run the dbg process.
4. Run the tip removal process.

### Preparing the sequencing reads

Let's prepare the sequencing reads.

`WorkSpaces` store sequencing reads in `ReadDatastore`s, provided by the
[ReadDatastores.jl](https://github.com/BioJulia/ReadDatastores.jl) package.

GenomeGraphs uses (and re-exports types and methods from) `ReadDatastores`.
`ReadDatastores` is a standalone package in its own right (although it was built
in the first place for GenomeGraphs).

If you want to use `ReadDatastores` as a standalone package in another Bio(Julia)
based project (and we recommend you do - the data stores are more efficient than
text files), you can find standalone docs for the package
[here](https://biojulia.net/ReadDatastores.jl/latest/). Some types and methods
documented there, are repeated in the Library section of this manual here, for
convenience.

Anyway, let's see how to build a paired end reads datastore!

```julia
using GenomeGraphs

fwq = open(FASTQ.Reader, "test/ecoli_pe_R1.fastq")
rvq = open(FASTQ.Reader, "test/ecoli_pe_R2.fastq")

ds = PairedReads(fwq, rvq, "ecoli-test-paired", "my-ecoli", 250, 300, 0, FwRv)
```

Here "ecoli-test-paired" is provided as the base filename of the datastore, the
datastore is given the name of "my-ecoli", this name will be used to identify it
in the workspace later.

The minimum length for the reads is set at 250 base
pairs, and the maximum length is set to 300 base pairs. Reads that are too short
will be discarded, reads that are too long are truncated.

!!! note
    The insert size of the paired reads to 0, since I'm not sure of it and
    right now the value is optional.

I set the orientation of the paired reads to `FwRv`. This is the default, and
means for every pair of reads, read 1 is oriented in the forward direction, and
read 2 is oriented backwards (forwards on the opposite strand). This orientation
distinguishes regular paired-end reads from other paired read types like
Long Mate Pairs.

### Add datastore to the WorkSpace

Now the datastore is created, it can be added to a workspace.

```julia
ws = WorkSpace()
add_paired_reads!(ws, ds)
```

### Run the `dbg` process

GenomeGraphs comes with some very high-level methods in it's API, that we like
to call *processes*. They perform some critical and common task as part of a
larger workflow. Examples include constructing a de-Bruijn graph from sequencing
reads, mapping reads to a graph, kmer counting and so on.

Once a workspace has an attached read datastore, you can run the dbg process to
produce a first de-Bruijn graph of the genome.

```julia
dbg!(BigDNAMer{61}, 10, ws, "my-ecoli")
```

!!! warning
    Some steps of this process, especially the Kmer counting steps, may take a
    long time for big inputs.
    No parallelism or batching to disk is used currently (although it is planned).
    This process should take just about a minute for E.coli paired end reads with
    a decent coverage.
    So, be warned, for big stuff, performance may suuuuuuucck in these early days!

### Run the `remove_tips` process

Now you have a raw compressed de-Bruijn assembly graph. You can start to use it
for analyses, and also try to improve its structure and resolve parts of the
graph that represent error, repetitive content, and so forth. Some of these
structures can be identified and resolved using only the topology of the graph,
some require additional information sources (linked reads / long reads / and so on),
to be incorporated into the workspace first.

Here let's see how to resolve a common structure, using only the graph topology.

Let's fix the tips of the graph.

Tips looks like this:

![tips](assets/tips.jpeg)

See how a piece of the assembly which should be one long stretch of sequence is
broken into 3 pieces (red) because of the existence of two tips (blue). Such
tips are defined topologically as very short segments which have one incoming
neighbour, and no outgoing neighbour.

Tips are caused by sequencer errors that occur at the end of reads, because the
sequencing by synthesis technique becomes more error prone over time; reagents
are consumed and products generated as time progresses, making the base detection
more difficult. Hence errors occur at the ends of reads, and erroneous kmers
from the read ends are unlikely to have forward neighbours, and they end up
forming tip nodes when they are incorporated into the de-Bruijn graph. 

you can remove these tips and improve the contiguity of the graph by using the
`remove_tips!` process.

```julia
remove_tips!(ws, 200)
```

The value of 200 provided is a size (in base pairs) threshold:
Nodes are considered tips only if they have one ingoing neighbour, no outgoing
neighbours, and are (in this case) smaller than 200 base pairs in length.

Once this process is finished you will have a collapsed de-Bruijn graph with the
tips removed.

## Using the NodeView interface

Graphs can be complicated data structures. If you want an in depth explanation
of the data-structure used to represent genome graphs, then head [here](notyet).

However, most people should not have to care about the internal structure of the
graph.

To make things as simple as possible, the `NodeView` type provides a single
entry point for node-centric analyses. The `NodeView` wraps a point to a
workspace's graph and contains a node id. A `NodeView` gives you acces to a
node's underlying sequence, the nodes neighbouring nodes in the forward and
backward directions, the reads mapped to a node, and kmer coverage over the node.

To get a `NodeView` of a node, use the `node` method on a `WorkSpace`, providing
a node id number. A positive ID denotes a view of the node traversing it in the
forward (canonical) direction. A negative ID denotes a view of the node,
traversing it in the reverse complement direction.   

```julia
julia> n = node(ws, 3)
A view of a graph node (node: 3, graph: sdg):
  AAAAAACCTCCGCAACCCCATGTTTTCACATAACTGTTG…GCCATGACCGGCTGGCTGTCAGGCTGTCACTGATAATCA


```