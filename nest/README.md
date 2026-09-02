# d-SEAMS cage counts as PLUMED collective variables

This folder is the PLUMED-NEST egg for `DSEAMS_CAGES`: the input file,
the module it loads and where the module comes from.

- `plumed.dat`: three actions, one each for ice nucleation in water, a
  growing front in an NaCl brine, and the filled cages of a methane
  hydrate. Every component is a count without derivatives, for `PRINT`,
  `COMMITTOR` and analysis rather than biasing.
- The module is `libdseams_plumed.so`, built from this repository
  against the d-SEAMS engine (`pixi run build`, or `meson setup build
  && meson compile -C build`). `LOAD FILE=` takes the path to that
  library; PLUMED 2.9 or later.
- The counts are the ones the d-SEAMS 2 paper measures against CHILL+,
  template matching and Steinhardt order parameters on the same
  frames; the reproducibility package runs these inputs through LAMMPS
  on a supercooled mW liquid with an ice seed and on a TIP4P/2005
  ice-brine interface.

Atom ranges in `plumed.dat` are for the example systems of the paper
(4096 mW molecules; 3000 waters and 60 ions; 2944 waters and 256
methanes) and want editing for another system: `ATOMS` is one oxygen
per molecule, `IONS` one atom per ion, `GUESTS` one atom per guest.
