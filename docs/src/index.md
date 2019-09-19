# GenomeGraphs

[![Latest release](https://img.shields.io/github/release/BioJulia/GenomeGraphs.svg)](https://github.com/BioJulia/GenomeGraphs/releases/latest)
[![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/BioJulia/GenomeGraphs.jl/blob/master/LICENSE) 
[![Stable documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://biojulia.github.io/GenomeGraphs/stable)
[![Pkg Status](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)
[![Chat](https://img.shields.io/gitter/room/BioJulia/GenomeGraphs.svg)](https://gitter.im/BioJulia/GenomeGraphs)

_A graph based genomics framework for the julia/BioJulia ecosystem_

## Introduction

GenomeGraphs is designed to do one thing - provide a framework that makes it
simple for a human to work with genome graphs from scripts or interactively.

Graphs are the core representation used by genome assemblers to represent
genome sequence models constructed from reads. At the time of writing it is fair
to say that until recently, their use has been limited to the internals of
genome assemblers, which are often treated as black boxes that output a series
of flattened sequences in FASTA format.

The use of graphs has increased in recent years thanks to the GFA file format
and developments in genome variation graphs and sequence to graph mappers.

However, a lack of inter operation between graph-based tools, and limited tools
for downstream graph-based analysis, contribute to a perceived complexity which
maintains linear sequences and FASTA files as the typical unit of genomic sequence
exchange.

Such flattening of graph representations within pipelines with multiple steps,
that use different types of sequencing in an iterative fashion, produces
ever-longer linear genome sequences through an information loss process. 
As a result, genome assembly projects are prone to error propagation and
difficult to reproduce and control.

To address these problems for the BioJulia ecosystem, GenomeGraphs provides a
flexible framework for building and integrating information over genome graphs.

## Framework overview

This package implements a `SequenceDistanceGraph` type that represents sequences
as nodes and the adjacency between sequences in links/edges. Rather than work
directly with this graph data structure, you interact with a `WorkSpace`. A
`WorkSpace` associates a `SequenceDistanceGraph` with raw sequencing reads, 
sequence to graph mappings, and k-mer counts. The `WorkSpace` and the API 
provides a working environment that enables you to project different kinds of 
information over a graph, and navigate and analyse each node of a sequence graph.

Within the `WorkSpace`, you will find `DataStore` types permit random access to
short, linked, and long read sequences stored on disk in `BioSequences.jl`'s
native bit encoding. Each datastore has an associated `Mapper` in the workspace
that contains the output from mapping said reads onto the graph. `KmerCounts`
allow you to compute k-mer coverage over the graph from sequencing data, enabling
coverage analysis. Additional `DistanceGraph`s define alternative topologies over
`SequenceDistanceGraph` nodes. They are typically used to represent longer range
linkage information from various sequencing technologies.

Finally, a `NodeView` abstraction provides a proxy to a node, with methods to
navigate a graph and access mapped data.

This framework is intended to allow the user to expore genome graphs interactively
and to create processing methods for assembly or downstream analyses.

## TLDR; Package features

- Efficient on-disk (buffered) data stores for sequencing reads.
- A simple graph data-structure for representing assembled genomes.
- TODO: Transparent mapping of sequencing reads onto graphs.
- TODO: Kmer counts and coverage projection over genome graph nodes.
- Workspaces binding a genome graph, mapped sequences, kmer counts, and annotation.
- De-novo genome assembly utilities:
  - de-Bruijn graph construction, with tip-clipping & bubble-popping.
