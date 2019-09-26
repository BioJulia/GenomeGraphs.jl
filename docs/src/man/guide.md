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

## Creating a WorkSpace

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

## Constructing a de-Bruijn graph

Let's see how to do option number 1, and construct a de-Bruijn graph from
raw sequencing reads. This can be achieved with a few simple steps:

1. Prepare the sequencing reads.

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
in the workspace later. The minimum length for the reads is set at 250 base
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

