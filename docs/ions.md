# Ions

Ions stay off the hydrogen-bond network. When `IONS` is supplied,
`DSEAMS_CAGES` builds the `K`-nearest-neighbour graphs and seeded ice cage
assignment on the oxygen atoms in `ATOMS` alone, then classifies each ion
against the water molecules inside its first solvation shell.

## First-shell classification

For each atom in `IONS`, the action gathers every oxygen atom within
`ION_CUTOFF` (default `3.5` Å) using minimum-image periodic distances and
checks its seeded cage label (hexagonal or double-diamond cage, including
ring completion when `COMPLETE` is set; independent of `CHILL+`):

| Component | Ion state |
| --- | --- |
| `nionice` | Non-empty first water shell in which every oxygen carries a cage label. |
| `nionfront` | Mixed first water shell (at least one cage-labelled oxygen and at least one unlabelled oxygen). |
| `nionliq` | No cage-labelled oxygen in the first water shell (including an empty shell). |

Every ion falls into exactly one class, so
`nionice + nionfront + nionliq` always equals the number of atoms in
`IONS`. Without `IONS`, all three components stay `0`.

These counts measure brine rejection at a growing ice front in electrolyte
simulations, such as TIP4P/2005 water with the Madrid 2019 ion parameters.

## Example

The deck below reads Ångström coordinates (`UNITS LENGTH=A` with
`LENGTH_SCALE=1.0`) for 3000 water oxygens and 60 ions:

```text
UNITS LENGTH=A
LOAD FILE=./build/libdseams_plumed.so

brine: DSEAMS_CAGES ...
  ATOMS=1-3000
  IONS=3001-3060
  ION_CUTOFF=3.5
  LENGTH_SCALE=1.0
  COMPLETE
...

PRINT ...
  ARG=brine.nice,brine.nmax,brine.nionice,brine.nionfront,brine.nionliq
  STRIDE=500
  FILE=BRINE
...
```

For trajectories in nanometres (PLUMED's default length unit), omit
`UNITS LENGTH=A` and `LENGTH_SCALE=1.0` so the default `LENGTH_SCALE=10.0`
converts positions to Ångströms before applying `ION_CUTOFF=3.5`.

On a 512-site cubic-diamond oxygen lattice with two water sites replaced by
substitutional ions (510 oxygens, 2 ions), every oxygen in each ion's first
shell belongs to the cage network:

```text
#! FIELDS time brine.nice brine.nmax brine.nionice brine.nionfront brine.nionliq
 0.000000      510      510        2        0        0
```
