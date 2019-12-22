```@meta
CurrentModule = GenomeGraphs.Graphs
```

# API: The Graphs submodule

The Graphs submodule contains the core graph data-structure used in
GenomeGraphs.jl: The SequenceDistanceGraph.

It also defines some of the common methods used in the rest of the higher-level
parts GenomeGraphs.jl:

- To query nodes/sequences in the graph.
- To find certain topological structures.
- To edit graph structure.
- To write/load SequenceDistanceGraph's to/from file.

Some methods, particularly ones that access the sequences stored in the graph,
and to edit topology, are marked as unsafe and are not exported deliberately.

All types and methods defined in this module are documented here.

!!! note
    This is a reference of an internal sub-module's API for developers and
    experienced users. First ask yourself if what you need isn't covered by
    the higher-level WorkSpace API. 

## Types

```@docs
SequenceDistanceGraph
SDG
SDGNode
SequenceDistanceGraphLink
```

## Public / Safe methods

### Graph nodes and sequences

```@docs
name
n_nodes
sequence
```

### Graph Topology

```@docs
source
destination
distance
is_forwards_from
is_backwards_from
find_link
forward_links
backward_links
get_next_nodes
get_previous_nodes
find_tip_nodes
find_tip_nodes!
find_all_unitigs
find_all_unitigs!
```

### Graph IO

```@docs
Graphs.write_to_gfa1
Graphs.load_from_gfa1!
```

## Internal / Unsafe methods

### Graph nodes and sequences

```@docs
Graphs.empty_seq
Graphs.is_deleted
Graphs.unsafe_sequence
Graphs.check_node_id
Graphs.nodes
Graphs.node_unsafe
Graphs.node
```

### Graph topology

```@docs
Graphs.links
Graphs.linksof_unsafe
```

### Graph editing and manipulation

```@docs
Graphs.add_node!
Graphs.remove_node!
Graphs.add_link!
Graphs.remove_link!
Graphs.disconnect_node!
Graphs.collapse_all_unitigs!
```
