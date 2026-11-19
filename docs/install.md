# Install

Build the shared module `libdseams_plumed.so` against an installed PLUMED
(2.9 or newer), or run the self-contained PLUMED-NEST source archive on
PLUMED 2.10 or newer.

Meson queries `plumed info --include-dir` for the PLUMED headers. The PLUMED
kernel resolves its symbols when `LOAD` opens the module, so
`libdseams_plumed.so` does not link against `libplumed`. An installed
`seams-core` (`pkg-config seams-core` reporting 2.14.0 or newer) supplies
the engine; otherwise Meson automatically downloads and builds the bundled
`seams-core` `v2.14.0` wrap.

## With pixi

[pixi](https://pixi.sh) installs PLUMED, compilers, Meson, Ninja, and the
engine dependencies from conda-forge:

```bash
git clone https://github.com/d-SEAMS/dseams-plumed
cd dseams-plumed
pixi install
pixi run build          # builds build/libdseams_plumed.so
pixi run test
```

`pixi run build` verifies the engine wrap pin, configures `build/`, and
compiles `build/libdseams_plumed.so`. `pixi run test` runs the four
`plumed driver` checks described in [Driver tests](tests.md).

## Without pixi

Put `plumed` (2.9 or newer) on `PATH` together with a C++20 compiler, Meson,
Ninja, `pkg-config`, Eigen 3.4, BLAS/LAPACK, Highway, and chemfiles. A Rust
toolchain is also required when falling back to the `seams-core` wrap (it is
not needed when a pre-installed `seams-core` 2.14.0 or newer is found via
`pkg-config`).

```bash
git clone https://github.com/d-SEAMS/dseams-plumed
cd dseams-plumed
meson setup build --buildtype=release -Dseams-core:with_tests=false
meson compile -C build
meson install -C build --prefix "$HOME/.local"
```

Passing `-Dseams-core:with_tests=false` keeps `meson test -C build` focused
on this module's driver tests.

## Load the module

`LOAD FILE=` takes the path to the compiled shared library:

```text
UNITS LENGTH=A
LOAD FILE=./build/libdseams_plumed.so

ice: DSEAMS_CAGES ATOMS=1-512 CUTOFF=3.5 LENGTH_SCALE=1.0 COMPLETE
PRINT ARG=ice.nice,ice.nmax STRIDE=1 FILE=ICE
```

Two loader rules matter in practice:

- **Path prefix**: `dlopen` does not search the current working directory
  for a bare filename. Always pass a relative path with a slash (such as
  `./build/libdseams_plumed.so`) or an absolute path.
- **Engine libraries**: when built via the Meson wrap,
  `build/libdseams_plumed.so` records a relative `RPATH` pointing into
  `build/subprojects/`. Point `LOAD FILE=` directly at
  `build/libdseams_plumed.so`, or install both the module and `seams-core`
  into the same prefix with `meson install -C build`. Copying only
  `build/libdseams_plumed.so` to another directory without installing its
  engine library will fail to load `libyodaLib.so`.

Check that PLUMED parses and loads the module from the repository root:

```bash
pixi run plumed driver \
  --natoms 512 --parse-only --kt 2.49 --plumed plumed.dat
```

## PLUMED-NEST egg

On PLUMED 2.10 or newer, the PLUMED-NEST egg
[plumID:26.018](https://www.plumed-nest.org/eggs/26/018/) (`v0.4.1`)
provides a self-contained source archive with `dseams_plumed.cpp`, the
required `seams-core` `v2.14.0` code, Eigen 3.4.0 headers, a ready-to-run
`plumed.dat`, `frame.xyz`, and reference outputs. No separate `seams-core`
build is needed: `LOAD FILE=dseams_plumed.cpp` invokes `plumed mklib` on
first use.

```bash
URL=https://github.com/d-SEAMS/dseams-plumed/releases/download/v0.4.1
curl -fLO "$URL/dseams-plumed-v0.4.1-nest.zip"
unzip dseams-plumed-v0.4.1-nest.zip
cd dseams-plumed-v0.4.1-nest
plumed driver --plumed plumed.dat --ixyz frame.xyz --length-units A
```

The command writes `ICE`, `BRINE`, and `HYDRATE`, matching `reference/ICE`,
`reference/BRINE`, and `reference/HYDRATE`.

!!! note "C++17 compiler requirement"
    `plumed mklib` must compile `dseams_plumed.cpp` with C++17. PLUMED 2.10
    configures C++17 by default; PLUMED 2.9 defaults to C++11 and compiles
    the archive only if PLUMED itself was built with C++17 enabled.
