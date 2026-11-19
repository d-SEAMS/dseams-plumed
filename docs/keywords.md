# Keyword reference

`DSEAMS_CAGES` evaluates topological cage, `CHILL+`, ion, topology-key, and
polyhedral-cage counts on a group of oxygen atoms at each step.

All distance cutoffs and radii are in Ångströms. `LENGTH_SCALE` converts
positions and periodic box vectors from the PLUMED length unit into
Ångströms (`10.0` for default nanometres, `1.0` after `UNITS LENGTH=A`).
Both orthorhombic and triclinic periodic boxes (in the LAMMPS convention
used by PLUMED) are supported.

## Keywords

| Keyword | Kind | Default | Meaning |
| --- | --- | --- | --- |
| `ATOMS` | atoms, required | | Oxygen atoms of the water molecules (at least one atom). |
| `CUTOFF` | compulsory | `3.5` | Neighbour cutoff in Å for the `CHILL+` graph. |
| `CANDIDATE` | compulsory | `5.5` | Candidate shell in Å for the `K`-nearest-neighbour search. |
| `K` | compulsory | `4` | Neighbours per molecule in the bonded graph. |
| `LENGTH_SCALE` | compulsory | `10.0` | Factor from the PLUMED length unit to Å (`10.0` for nm, `1.0` after `UNITS LENGTH=A`). |
| `COMPLETE` | flag | off | Fill the last vertex of six-rings whose other five vertices carry a cage label. |
| `IONS` | atoms | | Ions classified against the water cage assignment; excluded from the water graph. See [Ions](ions.md). |
| `ION_CUTOFF` | compulsory | `3.5` | Radius in Å of an ion's first water shell. |
| `HOPS` | compulsory | `3` | Bonds from the centre in each local topology key. |
| `LIBRARY` | optional | | One or more topology key library files, comma separated. See [Topology keys](topology.md). |
| `SIGNATURE` | optional | | Ring-size census (`4:6,6:8`) or named cage table. See [Cages and guests](cages.md). |
| `GUESTS` | atoms | | One atom per guest molecule placed in `SIGNATURE` cages; requires `SIGNATURE`. See [Cages and guests](cages.md). |
| `GUEST_RADIUS` | compulsory | `4.0` | Maximum distance in Å from a guest to its nearest `SIGNATURE` cage centre. |

Omit `IONS`, `LIBRARY`, `SIGNATURE`, and `GUESTS` when unused.

## Components

Every component is a non-periodic integer count with no derivatives. PLUMED
reserves underscores in component names, so names are written as single
lowercase words (for example `ice.chillmax`, not `ice.chill_max`).

### Ice cages and clusters

| Component | Meaning |
| --- | --- |
| `nice` | Molecules in a hexagonal or double-diamond cage (`nic + nih + nmixed`). |
| `nmax` | Largest connected cluster of cage molecules on the union graph. |
| `nclus` | Number of connected cage clusters on the union graph. |
| `nic` | Molecules in a double-diamond cage only (ice Ic). |
| `nih` | Molecules in a hexagonal cage only (ice Ih). |
| `nmixed` | Molecules belonging to both cage types. |

### CHILL+ and six-rings

| Component | Meaning |
| --- | --- |
| `chillice` | `CHILL+` bulk ice molecules (cubic plus hexagonal). |
| `chillmax` | Largest connected cluster of `CHILL+` bulk-ice molecules on the cutoff graph. |
| `chillinterfacial` | `CHILL+` interfacial ice molecules. |
| `sixrings` | Primitive six-membered rings on the union graph. |

### Ions

Without `IONS`, all three ion components stay `0`. When `IONS` is given,
`nionice + nionfront + nionliq` equals the number of atoms in `IONS`.

| Component | Meaning |
| --- | --- |
| `nionice` | Ions whose first water shell is non-empty and entirely cage-labelled. |
| `nionfront` | Ions with a mixed first water shell (some cage-labelled, some unlabelled). |
| `nionliq` | Ions with no cage-labelled water in the first shell (including empty shells). |

### Local topology keys

Without `LIBRARY`, `nnamed` stays `0`, while `nclasses` is always computed.

| Component | Meaning |
| --- | --- |
| `nclasses` | Distinct local topology classes on the mutual graph within `HOPS` bonds. |
| `nnamed` | Molecules whose local topology key appears in `LIBRARY`. |

### Signature cages and guests

Without `SIGNATURE`, `ncages` stays `0`. Without `GUESTS`, `noccupied`,
`nmultiple`, and `nfreeguest` stay `0`.

| Component | Meaning |
| --- | --- |
| `ncages` | Closed polyhedra matching `SIGNATURE`. |
| `noccupied` | `SIGNATURE` cages holding at least one guest (`noccupied <= ncages`). |
| `nmultiple` | `SIGNATURE` cages holding more than one guest (`nmultiple <= noccupied`). |
| `nfreeguest` | Guests farther than `GUEST_RADIUS` from every `SIGNATURE` cage centre. |

## Parse-time checks

`DSEAMS_CAGES` raises a PLUMED input error during setup if:

- `ATOMS` is empty.
- `GUESTS` is supplied without `SIGNATURE`.
- Any file listed in `LIBRARY` cannot be opened, or `LIBRARY` lists no file.
- The deepest library in `LIBRARY` was built at a hop count different from
  `HOPS`.
- Two files in `LIBRARY` share the same hop count.
- Files in `LIBRARY` disagree in colouring (`--color-types`).

## Minimal example

Count sodalite cages (`4:6,6:8`) on a periodic network of 324 oxygen sites:

```text
UNITS LENGTH=A
LOAD FILE=./build/libdseams_plumed.so

sod: DSEAMS_CAGES ...
  ATOMS=1-324
  CUTOFF=3.5
  LENGTH_SCALE=1.0
  SIGNATURE=sodalite
...

PRINT ARG=sod.ncages STRIDE=100 FILE=SOD
```
