# Topology keys

Every molecule carries a local topology key: its rooted bonded
neighbourhood on the mutual graph within `HOPS` bonds (three by default).
`nclasses` counts the distinct keys in the frame. Point `LIBRARY` at a
key library written by `seams fingerprint --emit-library LABEL`, or by
the Python and Lua front ends. Then `nnamed` counts the molecules whose
key the library names. The same count names a hydrate framework or a
high-pressure ice.

```text
hyd: DSEAMS_CAGES ATOMS=1-2944 LIBRARY=sI_sII.keys,sI_sII_h2.keys HOPS=3
PRINT ARG=hyd.nnamed,hyd.nclasses STRIDE=500 FILE=HYDRATE
```

`LIBRARY` takes several files, comma separated, built at different hop
counts. Each molecule takes the label of the deepest library that holds
its key, so a molecule whose wide neighbourhood is disturbed still gets
a name from its inner shells.

The deepest library must be built at `HOPS`. Two libraries at one hop
count, or libraries that differ in colouring, are errors. Without
`LIBRARY`, `nnamed` stays 0.
