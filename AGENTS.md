# Contributing to dseams-plumed

dseams-plumed builds `libdseams_plumed.so`, a PLUMED module that exposes
the d-SEAMS engine (seams-core) as actions: `DSEAMS_CAGES`,
`DSEAMS_CHILL`, `DSEAMS_RINGS` and `DSEAMS_SMOOTH_CHILL`.

## Layout

- `plumed/dseams/`: the action sources and `smooth_chill.hpp`, laid out
  as a PLUMED module directory (`module.type` default-off).
- `meson.build`: the module, the offline engine driver
  (`dseams_offline`), the C++ checks and the test registrations.
- `subprojects/seams-core.wrap`: the engine pin.
- `tests/`: driver checks in Python (`driver_check.py`,
  `chill_match.py`, `smooth_check.py`, `ion_check.py`, `guest_check.py`,
  `library_check.py`) and C++ checks (`thread_bits.cpp`,
  `track_check.cpp`, `offline_frame.cpp`).
- `regtest/`: PLUMED-style regtests. `regtest/run` copies each
  `regtest/dseams/rt-*` directory, runs `plumed driver`, and diffs every
  `*.reference` file.
- `examples/`: LAMMPS, GROMACS and OpenMM inputs, each with a README.
- `nest/`: the PLUMED-NEST egg.
- `docs/`: the action manual, fit and free-energy figures.
- `bench/logs/`: raw timing lines quoted in the README.
- `changelog.d/`: towncrier news fragments.

## Build and test

The pixi environment provides PLUMED, the compilers, meson and the
engine dependencies.

```bash
pixi run build   # meson setup + compile, build/libdseams_plumed.so
pixi run test    # meson test -C build --print-errorlogs
```

Without pixi, put `plumed` on `PATH` and run

```bash
meson setup build --buildtype=release --force-fallback-for=seams-core
meson compile -C build
meson test -C build --print-errorlogs
bash regtest/run   # DSEAMS_MODULE and PLUMED_BIN override the defaults
```

`driver_library` runs only when a `seams` binary is on `PATH`.

A change to an action keeps every regtest reference unchanged unless the
change is meant to move the number. Regenerate a `*.reference` only in
the commit that changes the behaviour, and say why in the message.
Never loosen a tolerance or skip a test to make a run pass.

## The seams-core pin

`subprojects/seams-core.wrap` names a seams-core revision. The
`project()` version in `meson.build` and the version in `pixi.toml`
follow the seams-core release the wrap pins: dseams-plumed 2.14.0 builds
against seams-core 2.14.0. The pixi `setup` task passes
`--force-fallback-for=seams-core`, so the wrap wins over an installed
engine.

To move the pin:

1. Point `revision` at the release tag's commit on
   `https://github.com/d-SEAMS/seams-core.git`.
2. Add `changelog.d/+engine-<rev>.changed.md` naming the revision.
3. Build and run the full test suite against the new engine.

## Commits and changelog

- Conventional Commits: `feat:`, `fix:`, `perf:`, `docs:`, `test:`,
  `chore:`, `refactor:`. One logical change per commit.
- Every user-visible change adds a fragment under `changelog.d/` in the
  same commit (`+slug.<type>.md`, types `added`, `changed`,
  `deprecated`, `removed`, `fixed`, `security`).
- The release cut is `towncrier build --version X.Y.Z`, which writes
  `CHANGELOG.md` and removes the consumed fragments.
- ASCII punctuation in sources, docs and messages.

## Data

Large outputs (trajectories, OPES `KERNELS`/`STATE` files, long
`COLVAR` tables) do not go in git. Anything over about 100 KB goes to a
dataset, and the repository keeps a manifest with the dataset path and
a checksum. The examples keep only the inputs and the reduced
free-energy tables their READMEs quote.
