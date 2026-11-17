# d-SEAMS cage counts as PLUMED collective variables

PLUMED 2.9 or later. The inputs were checked with PLUMED 2.9, PLUMED
2.10, and PLUMED master. The module loads with `LOAD FILE=` and needs
the d-SEAMS engine at the seams-core tag `v2.14.0`.

This folder is the PLUMED-NEST egg for `DSEAMS_CAGES`. `plumed.dat` is
the input. `nest.yml` is the submission metadata and is kept beside the
archive, not inside it.

- `plumed.dat`: three actions, one each for ice nucleation in water, a
  growing front in an NaCl brine, and the filled cages of a methane
  hydrate. Every component is a count without derivatives, for `PRINT`,
  `COMMITTOR` and analysis rather than biasing.
- The module is `libdseams_plumed.so`, built from
  https://github.com/d-SEAMS/dseams-plumed at tag `v0.4.0`
  (`pixi run build`, or `meson setup build && meson compile -C build`).
  The meson wrap fetches seams-core tag `v2.14.0`. An installed engine
  is accepted when `pkg-config seams-core` reports 2.14.0 or newer.
  `LOAD FILE=./libdseams_plumed.so` loads the copy in the working
  directory. A bare filename is not on the dynamic loader's search path.
- Authors: Rohit Goswami, Amrita Goswami, and Jayant K. Singh.
  Contributor: Rohit Goswami. Cite
  10.1021/acs.jcim.0c00031. The d-SEAMS 2.0 manuscript does not yet have
  a journal DOI.
- The counts are the ones the d-SEAMS 2 paper measures against CHILL+,
  template matching and Steinhardt order parameters on the same
  frames; the reproducibility package runs these inputs through LAMMPS
  on a supercooled mW liquid with an ice seed and on a TIP4P/2005
  ice-brine interface.

Atom ranges in `plumed.dat` are for the example systems of the paper
(4096 mW molecules; 3000 waters and 60 ions; 2944 waters and 256
methanes) and want editing for another system: `ATOMS` is one oxygen
per molecule, `IONS` one atom per ion, `GUESTS` one atom per guest.

NEST checks syntax with:

```
plumed driver --natoms 100000 --parse-only --kt 2.49 --plumed plumed.dat
```

That command dlopens `libdseams_plumed.so` from the working directory.
The archive does not vendor the library: it links the seams-core
engine, which the NEST runners do not have installed.
