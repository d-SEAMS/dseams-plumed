# dseams-plumed

`dseams-plumed` exposes the d-SEAMS topological ice, cage, ion, and local
environment classifiers as a PLUMED collective variable action named
`DSEAMS_CAGES`.

On each step, `DSEAMS_CAGES` reads a group of oxygen atoms and builds the
mutual and union four-nearest-neighbour graphs, the primitive six-membered
rings, and the seeded cage assignment of d-SEAMS 2.0. With `COMPLETE` it
fills the last vertex of six-rings whose other vertices carry a label.
Alongside cage occupancy, largest connected cluster, and the ice Ic / ice Ih
/ mixed split, the action reports the `CHILL+` classification on the same
frame, optional ion rejection counts, local topology key matches, and closed
polyhedral cage counts with guest occupancy.

!!! note "Counts without derivatives"
    Every component registered by `DSEAMS_CAGES` is an integer count with no
    derivatives. Use them with `PRINT`, `COMMITTOR` basins, stopping
    conditions, and `plumed driver` analysis, not with biasing actions such
    as `RESTRAINT` or `METAD`.

## Quick start

Lengths inside the engine are in Ångströms. `LENGTH_SCALE` converts from the
PLUMED length unit (`10.0` for the default nanometres, `1.0` after
`UNITS LENGTH=A`).

```text
UNITS LENGTH=A
LOAD FILE=./build/libdseams_plumed.so

ice: DSEAMS_CAGES ...
  ATOMS=1-4096
  CUTOFF=3.5
  LENGTH_SCALE=1.0
  COMPLETE
...

PRINT ...
  ARG=ice.nice,ice.nmax,ice.chillice,ice.chillmax
  STRIDE=100
  FILE=ICE
...

COMMITTOR ...
  ARG=ice.nmax
  STRIDE=100
  BASIN_LL1=0 BASIN_UL1=20
  BASIN_LL2=800 BASIN_UL2=100000
...
```

Omit `UNITS LENGTH=A` and `LENGTH_SCALE=1.0` when reading trajectories or
running engines that pass coordinates in nanometres.

## Where to go next

- [Install](install.md): build `libdseams_plumed.so` with pixi or Meson, or
  run the self-contained PLUMED-NEST egg
  [plumID:26.018](https://www.plumed-nest.org/eggs/26/018/) (`v0.4.1`).
- [Keyword reference](keywords.md): every keyword, default, component, and
  parse-time check in `DSEAMS_CAGES`.
- [Cages and guests](cages.md): count closed polyhedra by ring-size
  `SIGNATURE` (clathrate hydrates, sodalite, custom cages) and assign
  `GUESTS` to cage centres.
- [Ions](ions.md): classify ions into ice, liquid, and front shells for
  brine rejection studies.
- [Topology keys](topology.md): match local bonded neighbourhoods against
  single-hop or multi-hop `seams fingerprint` libraries.
- [Inside LAMMPS](lammps.md): `fix plumed` ordering, length units, and
  multi-site water atom selection.
- [Driver tests](tests.md): verify a build and run `plumed driver` on
  trajectories.

## Companion documentation and video

The main [d-SEAMS documentation](https://docs.dseams.info) covers the
companion `SEAMS` PLUMED action built inside `seams-core` (which provides a
differentiable `q6` component alongside cubic and ring fractions) and links
to this module:

- [PLUMED how-to](https://docs.dseams.info/howto/engine-plumed.html)
- [PLUMED walkthrough](https://docs.dseams.info/walkthrough/plumed.html)

The walkthrough video is available on
[Odysee](https://odysee.com/dseams-plumed:ddd2f9e27e0a126ef3feba3bb6da8378464e57df).

<div class="video" markdown="0">
<iframe
  src="https://odysee.com/$/embed/dseams-plumed:ddd2f9e27e0a126ef3feba3bb6da8378464e57df"
  title="dseams-plumed walkthrough"
  width="100%" style="aspect-ratio: 16 / 9;"
  loading="lazy" allowfullscreen></iframe>
</div>

## Related projects

- [seams-core](https://github.com/d-SEAMS/seams-core): the C++ engine and
  the `seams` CLI
- [PydSEAMSlib](https://github.com/d-SEAMS/PydSEAMSlib): Python front end
- [yodaStruct](https://github.com/d-SEAMS/yodaStruct): Lua and Fennel front end
- [dseams.info](https://dseams.info): project home page
