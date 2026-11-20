# Driver tests

Meson registers four `plumed driver` checks that exercise every mode of
`DSEAMS_CAGES`, plus `feature_query`. When the build uses the seams-core
wrap with both structure readers, it also registers four `seams read`
failure checks. Run the full suite after building:

```bash
pixi run test
```

For a non-pixi build, run Meson directly:

```bash
meson test -C build --print-errorlogs
```

## Registered checks

### Nucleation (`driver_nucleation`)

`tests/driver_check.py` converts a LAMMPS dump trajectory into an xyz file
with per-frame periodic cell vectors on the comment line, loads
`libdseams_plumed.so`, and runs `plumed driver` without requiring the
molfile plugin.

With no trajectory arguments, the script tests the bundled 4096-molecule mW
cubic ice fixture (`mW_cubic.lammpstrj`) and verifies that every molecule
belongs to a single cage cluster (`nice == nmax == 4096`):

```bash
pixi run python tests/driver_check.py \
  "$(pixi run which plumed)" "$PWD/build/libdseams_plumed.so"
```

You can also validate against your own LAMMPS dump and a `walk_compare`
reference table from `seams-core` to check `nice`, `nmax`, `chillice`, and
`chillmax` frame by frame:

```bash
pixi run python tests/driver_check.py \
  "$(pixi run which plumed)" "$PWD/build/libdseams_plumed.so" \
  trajectory.lammpstrj walk_compare.txt
```

### Ions (`driver_ions`)

`tests/ion_check.py` builds a 512-site cubic-diamond oxygen lattice with two
sites replaced by substitutional ions (510 oxygens, 2 ions) and checks
`nionice = 2`, `nionfront = 0`, `nionliq = 0`, while the remaining oxygen
lattice retains its cage assignment.

### Guests (`driver_guests`)

`tests/guest_check.py` builds a periodic 5x5x5 simple-cubic lattice (125
sites) with `K=6`, `SIGNATURE=4:6`, and 64 guests at unit-cube centres, and
checks `ncages = 125`, `noccupied = 64`, `nmultiple = 0`, and
`nfreeguest = 0`.

### Library (`driver_library`)

`tests/library_check.py` writes a cubic-diamond lattice, emits 3-hop and
2-hop libraries with `seams fingerprint --emit-library`, runs `plumed driver`
with `LIBRARY=ic.keys,ic2.keys`, and checks `nclasses = 1` and
`nnamed = 512`. Meson locates `seams` on `PATH` or uses the CLI built by the
`seams-core` wrap.

### Feature query (`feature_query`)

`DSEAMS_CAGES` never opens a structure file (frames come from the MD
engine). `tests/feature_check.py` checks the two `extern "C"`
functions the built `libdseams_plumed.so` exports, reporting whether the
linked `seams-core` engine was built with its optional
[chemfiles](https://chemfiles.org) (PDB/GRO/DCD) and
[readcon-core](https://github.com/lode-org/readcon-core) (`.con`) readers:

```bash
pixi run python tests/feature_check.py \
  "$(pixi run which plumed)" "$PWD/build/libdseams_plumed.so"
```

The script preloads the PLUMED kernel with `RTLD_GLOBAL` before loading the
module, the same resolution `LOAD` gets inside a running `plumed driver`,
then calls `dseams_plumed_has_chemfiles()` and `dseams_plumed_has_readcon()`
through `ctypes`. With `--require-chemfiles`/`--require-readcon` the check
fails (does not skip) when that reader is absent; CI's `build-test` job
configures with `-Dseams-core:with_chemfiles=enabled
-Dseams-core:with_readcon=enabled` explicitly (a deliberate configuration,
not relying on auto-detection), so it runs the script that way to make
sure the assertion cannot pass silently.

### Missing and unreadable structure files (`reader_fail_pdb`, `reader_fail_gro`, `reader_fail_dcd`, `reader_fail_con`)

Registered only when this build uses the seams-core wrap and both readers
are compiled in. `tests/reader_fail_check.py` runs `seams read` on a
nonexistent file, an empty file, and a mode-000 file for one extension
(`pdb`, `gro`, `dcd`, or `con`). Each run must exit 2, write nothing to
stdout, and print a stderr line `cannot read <path>`. The tests do not
read a file that contains atoms.

```bash
pixi run python tests/reader_fail_check.py \
  "$PWD/build/subprojects/seams-core/seams" pdb
```

## Running `plumed driver` on your own xyz trajectory

PLUMED's `--ixyz` reader accepts nine numbers on each frame's comment line
for the three lattice vectors (`ax ay az bx by bz cx cy cz`), allowing
per-frame NPT box changes without `--box`. Pair `--length-units A` with
`UNITS LENGTH=A` and `LENGTH_SCALE=1.0` in `plumed.dat`:

```bash
plumed driver --plumed plumed.dat --ixyz traj.xyz --length-units A
```

To verify syntax and module loading without a trajectory:

```bash
plumed driver --natoms 4096 --parse-only --kt 2.49 --plumed plumed.dat
```
