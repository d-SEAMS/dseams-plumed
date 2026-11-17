# Driver tests

Meson registers four `plumed driver` checks; `pixi run test` runs them.

## Nucleation

`tests/driver_check.py` converts a LAMMPS dump to xyz and runs the
driver with the module loaded. Per frame it compares `nice`, `nmax`,
`chillice` and `chillmax` with the engine `walk_compare` table. Absent
a dump, the script uses the mW cubic fixture under `DSEAMS_PLUMED_DATA`
or the seams-core subproject. The dump becomes xyz with the cell on the
comment line. The check does not depend on the molfile plugin.

```bash
python tests/driver_check.py $(which plumed) build/libdseams_plumed.so \
    data/nucleation.lammpstrj results/walks/figshare-nucleation.txt
```

## Ions

`driver_ions` (`tests/ion_check.py`) places two ions on a cubic-diamond
oxygen lattice. It checks `nionice` 2, `nionfront` 0, `nionliq` 0, with
`nice` still the remaining lattice.

## Guests

`driver_guests` (`tests/guest_check.py`) builds a 5x5x5 simple-cubic
lattice with `SIGNATURE=4:6` and a guest at 64 of the 125 cube centres.
It checks `ncages` 125, `noccupied` 64, `nmultiple` 0 and
`nfreeguest` 0.

## Library

`driver_library` (`tests/library_check.py`) writes a cubic-diamond
lattice, emits its key library with `seams fingerprint --emit-library`,
and checks that `nnamed` equals the molecule count and `nclasses` is
one. Meson takes `seams` from `PATH`, or the CLI the seams-core wrap
builds.

## Continuous integration

Every push to `main` checks that the wrap revision is the seams-core
tag `v2.14.0`, builds the module in the pixi environment, and runs the
driver tests.
