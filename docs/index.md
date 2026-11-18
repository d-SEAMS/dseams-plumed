# dseams-plumed

The d-SEAMS topological ice score as a PLUMED action.

`DSEAMS_CAGES` takes a group of oxygen atoms and, every step, builds the
mutual and union four-nearest-neighbour graphs, the primitive six-membered
rings and the seeded cage assignment of d-SEAMS 2.0. With `COMPLETE` it
fills the last vertex of six-rings whose other vertices carry a label.
The action reports cage occupancy, the largest connected cluster, the
split into ice Ic, ice Ih and mixed, and the four-neighbour star method
`CHILL+` on the same frame. Integer counts without derivatives drive
`PRINT`, `COMMITTOR` basins and stopping conditions. Cage lengths
inside the engine use the unit A. `LENGTH_SCALE` converts from the
PLUMED length unit (10 for the default nm, 1 after `UNITS LENGTH=A`).

```text
LOAD FILE=./libdseams_plumed.so
ice: DSEAMS_CAGES ATOMS=1-4096 CUTOFF=3.5 COMPLETE
PRINT ARG=ice.nice,ice.nmax,ice.chillice,ice.chillmax STRIDE=100 FILE=ICE
COMMITTOR ARG=ice.nmax STRIDE=100 BASIN_LL1=0 BASIN_UL1=20 BASIN_LL2=800 BASIN_UL2=100000
```

[Install](install.md) builds `libdseams_plumed.so`, or use the
PLUMED-NEST egg [plumID:26.018](https://www.plumed-nest.org/eggs/26/018/)
(v0.4.1). The [keyword reference](keywords.md) lists every keyword and
component the action registers.

## Book and video

The [d-SEAMS documentation](https://docs.dseams.info) covers the module:

- [PLUMED how-to](https://docs.dseams.info/howto/engine-plumed.html)
- [PLUMED walkthrough](https://docs.dseams.info/walkthrough/plumed.html)

The walkthrough video is on
[Odysee](https://odysee.com/dseams-plumed:ddd2f9e27e0a126ef3feba3bb6da8378464e57df).

<div class="video" markdown="0">
<iframe
  src="https://odysee.com/$/embed/dseams-plumed:ddd2f9e27e0a126ef3feba3bb6da8378464e57df"
  title="dseams-plumed walkthrough"
  width="100%" style="aspect-ratio: 16 / 9; border: 0;"
  loading="lazy" allowfullscreen></iframe>
</div>

## Related projects

- [seams-core](https://github.com/d-SEAMS/seams-core): the C++ engine and
  the `seams` CLI
- [PydSEAMSlib](https://github.com/d-SEAMS/PydSEAMSlib): Python front end
- [yodaStruct](https://github.com/d-SEAMS/yodaStruct): Lua and Fennel front end
- [dseams.info](https://dseams.info): project site
