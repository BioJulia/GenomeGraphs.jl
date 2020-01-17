```@meta
CurrentModule = GenomeGraphs.GraphIndexes
```

# API: The GraphIndexes submodule

The GraphIndexes submodule contains some simple indexes and data structures that
are useful if not required, for the development of mappers and motif/sequence
finding heuristics.

All types and methods defined in this module are documented here.

!!! note
    This is a reference of an internal sub-module's API for developers and
    experienced users. First ask yourself if what you need isn't covered by
    the higher-level WorkSpace API.

!!! warning
    This submodule is still under construction.

## Types

```@docs
UniqueMerIndex
```

## Public / Safe methods

```@docs
isunmappable
find_unique_mer
```

## Internal / Unsafe methods