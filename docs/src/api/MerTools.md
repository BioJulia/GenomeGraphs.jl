```@meta
CurrentModule = GenomeGraphs.MerTools
```

# API: The MerTools submodule



!!! note
    This is a reference of an internal sub-module's API for developers and
    experienced users. First ask yourself if what you need isn't covered by
    the higher-level WorkSpace API. 

## Types

```@docs
MerCount
MerCountHist
DNAMerCount
RNAMerCount
```

## Public / Safe methods

```@docs
mer
freq
collapse_into_counts
collapse_into_counts!
merge_into!
build_freq_list
```

## Internal / Unsafe methods

```@docs
unsafe_collapse_into_counts!
unsafe_merge_into!
```