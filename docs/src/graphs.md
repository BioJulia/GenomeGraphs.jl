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

```@docs
add_node!
remove_node!
add_link!
remove_link!
disconnect_node!
```