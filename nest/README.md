# d-SEAMS cage counts as PLUMED collective variables

This archive is self-contained source for PLUMED 2.9 or later.
`LOAD FILE=dseams_plumed.cpp` asks PLUMED to compile the action on first
use. The generated source contains the required seams-core v2.14.0 code.
The `Eigen` directory supplies its header-only matrix types; no d-SEAMS
library is required.

`plumed.dat` uses only the core `UNITS`, `LOAD`, and `PRINT` directives.
`DSEAMS_CAGES` is defined by the bundled source. It does not require the
optional symmetry-functions, environment-similarity, or adjacency-matrix
modules.

The compiler used by `plumed mklib` must support C++17. PLUMED 2.10
already configures that language level. PLUMED 2.9 builds must likewise
be configured with C++17.

The example contains a cubic-diamond ice lattice, a second lattice with
two substitutional ions, and a periodic cubic cage lattice with 64
guests. The three actions write `ICE`, `BRINE`, and `HYDRATE`.

Run the complete example from this directory:

```sh
plumed driver --plumed plumed.dat --ixyz frame.xyz --length-units A
```

Compare the data rows with `reference/ICE`, `reference/BRINE`, and
`reference/HYDRATE`. PLUMED may also leave a versioned shared object
beside the source; it is a local compilation cache and is not needed in
the archive.

The PLUMED-NEST syntax check is:

```sh
plumed driver --natoms 100000 --parse-only --kt 2.49 --plumed plumed.dat
```

`ATOMS` is one oxygen per water, `IONS` is one atom per ion, and
`GUESTS` is one atom per guest. Edit the ranges and parameters in
`plumed.dat` for another system.

Authors: Rohit Goswami, Amrita Goswami, and Jayant K. Singh. Cite
10.1021/acs.jcim.0c00031.
