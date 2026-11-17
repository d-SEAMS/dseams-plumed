# Keyword reference

`DSEAMS_CAGES` registers the keywords and components below
in `registerKeywords` of `src/DseamsCages.cpp`.
Lengths are in A; `LENGTH_SCALE` converts the PLUMED length unit.

## Keywords

| Keyword | Kind | Default | Meaning |
| --- | --- | --- | --- |
| `ATOMS` | atoms, required | | The oxygen atoms of the water molecules. |
| `CUTOFF` | compulsory | `3.5` | Neighbour cutoff in A for the `CHILL+` graph. |
| `CANDIDATE` | compulsory | `5.5` | Candidate shell in A for the four-nearest walk. |
| `K` | compulsory | `4` | Neighbours per molecule in the bonded graph. |
| `LENGTH_SCALE` | compulsory | `10.0` | Factor from the PLUMED length unit to A (10 for nm, 1 after `UNITS LENGTH=A`). |
| `COMPLETE` | flag | off | Fill the last vertex of six-rings whose other vertices carry a label. |
| `IONS` | atoms | | Ions read against the water assignment; not part of the graph. See [Ions](ions.md). |
| `ION_CUTOFF` | compulsory | `3.5` | Radius in A of an ion's first water shell. |
| `HOPS` | compulsory | `3` | Bonds from the centre in each local topology key. |
| `LIBRARY` | optional | | One or more topology key libraries, comma separated. See [Topology keys](topology.md). |
| `SIGNATURE` | optional | | A ring-size census (`4:6,6:8`) or a named table. |
| `GUESTS` | atoms | | One atom per guest molecule; needs `SIGNATURE`. |
| `GUEST_RADIUS` | compulsory | `4.0` | A guest belongs to the nearest cage centre within this distance (A). |

`SIGNATURE` takes the named tables `sodalite`, `alpha`, `512`, `51262`,
`51264`, `51268`, `sh`, `sH`, `hc` and `ddc`. `hc` and `ddc` keep the
seeded ice finders. `ncages` counts the closed polyhedra that match.

`GUESTS` lists one atom per guest molecule (methane carbon, THF oxygen,
an ion). Each guest goes to the nearest cage centre, the periodic
centroid of a cage's vertices, within `GUEST_RADIUS`. `GUESTS` without
`SIGNATURE` is an error.

`LIBRARY` takes files written by `seams fingerprint --emit-library`.
Several files, built at different hop counts, name each molecule by the
deepest library that holds its key. The deepest library must be built at
`HOPS`. Two libraries at one hop count, or libraries that differ in
colouring, are errors.

```text
sod: DSEAMS_CAGES ATOMS=1-324 CUTOFF=3.5 SIGNATURE=sodalite
PRINT ARG=sod.ncages STRIDE=100 FILE=SOD
```

Omit `LIBRARY`, `IONS`, `SIGNATURE` or `GUESTS` when unused.

## Components

Each component is an integer count with no derivative. Use them in
`PRINT`, `COMMITTOR` and analysis, not in a bias.

| Component | Count |
| --- | --- |
| `nice` | Molecules in a hexagonal or double-diamond cage. |
| `nmax` | Largest connected cluster of cage molecules. |
| `nclus` | Number of connected cage clusters. |
| `nic` | Molecules in a double-diamond cage only. |
| `nih` | Molecules in a hexagonal cage only. |
| `nmixed` | Molecules in both cage types. |
| `chillice` | `CHILL+` cubic plus hexagonal molecules. |
| `chillmax` | Largest `CHILL+` bulk-ice cluster. |
| `chillinterfacial` | `CHILL+` interfacial molecules. |
| `sixrings` | Primitive six-membered rings on the union graph. |
| `nionice` | Ions whose first water shell is all cage molecules. |
| `nionfront` | Ions with a mixed first shell. |
| `nionliq` | Ions with no cage molecule in the first shell. |
| `nclasses` | Distinct local topology classes on the mutual graph. |
| `nnamed` | Molecules whose local key the `LIBRARY` names. |
| `ncages` | Closed polyhedra matching `SIGNATURE`. |
| `noccupied` | `SIGNATURE` cages holding at least one guest. |
| `nmultiple` | `SIGNATURE` cages holding more than one guest. |
| `nfreeguest` | Guests in no `SIGNATURE` cage. |

Without `IONS` the three ion components stay 0. Without `LIBRARY`,
`nnamed` stays 0. Without `SIGNATURE`, `ncages` stays 0; without
`GUESTS`, the three guest components stay 0.
