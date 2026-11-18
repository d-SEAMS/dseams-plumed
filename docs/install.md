# Install

Build the shared module against an installed PLUMED (2.9 or newer).
Meson runs `plumed info --include-dir` for the headers. The kernel
resolves the symbols when `LOAD` dlopens the file, so the module does
not link against libplumed. An installed `seams-core`
(`pkg-config seams-core`) supplies the engine; the wrap is the fallback.

The wrap pins the seams-core tag `v2.14.0`. An installed engine is
accepted when `pkg-config seams-core` reports 2.14.0 or newer. The
setup passes `-Dseams-core:with_tests=false` so `meson test` runs this
module's driver tests. The engine's own suite stays in seams-core.

## With pixi

[pixi](https://pixi.sh) supplies PLUMED, the compilers and meson from
conda-forge:

```bash
git clone https://github.com/d-SEAMS/dseams-plumed
cd dseams-plumed
pixi install
pixi run build          # build/libdseams_plumed.so
pixi run test
```

`pixi run build` checks the engine pin, configures `build/` and compiles.
`pixi run test` runs the driver tests described under
[Driver tests](tests.md).

## Without pixi

Put an installed `plumed` on `PATH`, with a C++ compiler, meson, ninja,
pkg-config, Eigen 3.4, BLAS/LAPACK, Highway, chemfiles and a Rust
toolchain for the engine, and run:

```bash
meson setup build --buildtype=release -Dseams-core:with_tests=false
meson compile -C build
meson install -C build   # libdseams_plumed.so into libdir
```

## Load the module

`LOAD FILE=` takes the path of the built library:

```text
LOAD FILE=./libdseams_plumed.so
```

A bare filename is not searched in the working directory; use
`./libdseams_plumed.so` or an absolute path.

## PLUMED-NEST egg

On PLUMED 2.10 or newer, the PLUMED-NEST egg
[plumID:26.018](https://www.plumed-nest.org/eggs/26/018/) (v0.4.1)
bundles a self-contained `dseams_plumed.cpp` with the required
`seams-core` `v2.14.0` code and Eigen 3.4.0 headers. Without a separate
`seams-core` build, `LOAD FILE=dseams_plumed.cpp` compiles the action on
first use, and the archive includes `plumed.dat` and `frame.xyz`.
