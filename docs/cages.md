# Cages and guests

`SIGNATURE` enumerates closed polyhedra on the union `K`-nearest-neighbour
graph by their ring-size face census. Adding `GUESTS` assigns guest
molecules (such as methane, THF, or ions) to those cages so you can track
filled-cage counts and multiple occupancy directly inside PLUMED.

## Cage signatures

Pass either a custom ring-size census of the form `size:count,...` (for
example `4:6` for a cube or `4:6,6:8` for a truncated octahedron) or one of
the built-in named tables. Named tables are matched case-insensitively.

| Name | Ring-size census | Polyhedron / role |
| --- | --- | --- |
| `sodalite` | `4:6,6:8` | Truncated octahedron (sodalite cage) |
| `alpha` | `4:6,6:8,8:6` | Alpha cage (LTA) |
| `512` | `5:12` | Pentagonal dodecahedron (5<sup>12</sup>, small sI / sII / sH cage) |
| `51262` | `5:12,6:2` | Tetrakaidecahedron (5<sup>12</sup>6<sup>2</sup>, large sI cage) |
| `51264` | `5:12,6:4` | Hexakaidecahedron (5<sup>12</sup>6<sup>4</sup>, large sII cage) |
| `51268` | `5:12,6:8` | Icosahedral 5<sup>12</sup>6<sup>8</sup> cage (large sH cage) |
| `sh`, `sH` | `4:3,5:6,6:3` | Medium sH cage (4<sup>3</sup>5<sup>6</sup>6<sup>3</sup>) |
| `hc` | `4:6,6:2` | Hexagonal ice cage (uses the seeded ice finder) |
| `ddc` | `6:7` | Double-diamond ice cage (uses the seeded ice finder) |

`hc` and `ddc` keep the seeded ice cage finders, while all other signatures
enumerate closed polyhedra on the primitive ring network up to the largest
ring size in `SIGNATURE`.

Match `K` and `CANDIDATE` to the host lattice coordination:

- Tetrahedral water and clathrate hydrate networks use the defaults `K=4`
  and `CANDIDATE=5.5` Å.
- Six-coordinated simple-cubic networks (for example `SIGNATURE=4:6`) use
  `K=6` and a `CANDIDATE` shell slightly above the lattice spacing.

## Guest occupancy

`GUESTS` takes one atom per guest molecule (for example, methane carbon or
THF oxygen). `GUESTS` requires `SIGNATURE`; omitting `SIGNATURE` raises a
setup error.

For each cage matching `SIGNATURE`, the action computes the cage centre as
the periodic centroid of its vertices (unwrapping every vertex to its
minimum image relative to the first vertex). Each guest atom is assigned to
the nearest cage centre within `GUEST_RADIUS` (default `4.0` Å). Guests
farther than `GUEST_RADIUS` from every cage centre count in `nfreeguest`.

| Component | Meaning |
| --- | --- |
| `ncages` | Closed polyhedra matching `SIGNATURE`. |
| `noccupied` | Matching cages containing at least one guest. |
| `nmultiple` | Matching cages containing two or more guests. |
| `nfreeguest` | Guests not within `GUEST_RADIUS` of any matching cage centre. |

## Clathrate hydrate example

Track small 5<sup>12</sup> cages and methane carbon occupancy in a hydrate
frame with Ångström coordinates:

```text
UNITS LENGTH=A
LOAD FILE=./build/libdseams_plumed.so

hyd: DSEAMS_CAGES ...
  ATOMS=1-2944
  GUESTS=2945-3456
  CUTOFF=3.5
  CANDIDATE=5.5
  K=4
  LENGTH_SCALE=1.0
  SIGNATURE=512
  GUEST_RADIUS=4.0
...

PRINT ...
  ARG=hyd.ncages,hyd.noccupied,hyd.nmultiple,hyd.nfreeguest
  STRIDE=500
  FILE=HYDRATE
...
```

## Simple-cubic cage example

On a periodic 5x5x5 simple-cubic lattice (125 host atoms, spacing 3.0 Å)
with 64 guests placed at unit-cube centres, every unit cube is a closed
`4:6` polyhedron:

```text
UNITS LENGTH=A
LOAD FILE=./build/libdseams_plumed.so

cube: DSEAMS_CAGES ...
  ATOMS=1-125
  GUESTS=126-189
  CUTOFF=3.5
  CANDIDATE=4.5
  K=6
  LENGTH_SCALE=1.0
  SIGNATURE=4:6
  GUEST_RADIUS=3.0
...

PRINT ...
  ARG=cube.ncages,cube.noccupied,cube.nmultiple,cube.nfreeguest
  STRIDE=1
  FILE=CUBES
  FMT=%8.0f
...
```

Running this deck on the 5x5x5 lattice yields:

```text
#! FIELDS time cube.ncages cube.noccupied cube.nmultiple cube.nfreeguest
 0.000000      125       64        0        0
```
