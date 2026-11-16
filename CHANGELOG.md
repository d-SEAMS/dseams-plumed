# Changelog

All notable changes to this project are documented in this file.

The format is [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Unreleased notes live in [`changelog.d/`](changelog.d/) and are assembled
by [towncrier](https://towncrier.readthedocs.io/).

<!-- towncrier release notes start -->

## [2.14.0] - 2026-11-15

### Added

- An OpenMM TIP4P-Ew example loads the module through PlumedForce and
  runs a short periodic box.
- DSEAMS_CAGES reports noccupiedinside, nmultipleinside and
  nfreeguestinside, the guest counts from the engine's ray-parity inside
  test, beside the centroid-within-radius counts.
- DSEAMS_CHILL reports the published Nguyen and Molinero 2015 bulk-ice
  fraction, and RULE=ovito-chill+ selects the OVITO 3.17.2 modifier. The
  labels have no derivatives. DSEAMS_RINGS reports the primitive ring
  census on the cutoff graph, also without derivatives.
  DSEAMS_SMOOTH_CHILL is a differentiable surrogate of that fraction for
  metadynamics and other biases. Its neighbour list is reused across
  steps, and the Y_3m gradient is analytic, so one evaluation on 4096
  molecules costs about the same as the discrete count.
- CITATION.cff with the d-SEAMS and CHILL+ papers.
- pixi environments `plumed29` and `plumed210`, and a committed
  `pixi.lock`. CI runs the full suite in both, on `main`, `release/**`
  and `v*` tags, and checks that the versions agree.

### Changed

- DSEAMS_SMOOTH_CHILL evaluates on PLUMED's NUM_THREADS. The pair list
  order is unchanged, so one thread and several threads return the same
  value, forces and virial.
- Engine wrap follows the seams-core 2.14.0 release (628b0af9), from
  v2.9.2 by way of v2.10.0 and the main revisions after it.
- The module is built for the toolchain's generic CPU, so it loads on
  any machine the compiler targets; `-Dnative=true` adds
  `-march=native`.
- The engine's own LAMMPS, GROMACS, OpenMM and PLUMED hooks are not
  built under the module, so a LAMMPS in the environment no longer
  breaks the build.
- The GROMACS slab example runs eight walkers on one shared OPES bias
  with the capped smooth fraction, and `scripts/fes_blocks.py` gives the
  block error of the free-energy difference.
- The TIP4P/Ice slab free energy is reported at 0.5, 1 and 2 ns, with
  block standard errors. At 330 K the minimum sits at low ice fraction.
- The plugin version follows the seams-core release it wraps, so this
  release is 2.14.0.
- The smooth fraction gathers Y_3m from a molecule-major neighbour list.
  On 32768 molecules the steady step is 5.3 times faster on 8 threads
  than on 1, and the value, forces and virial stay bit-identical.
- The smooth ice fraction follows the published CHILL+ counts: four
  staggered neighbours, or three staggered and one eclipsed. On a ladder
  from liquid to ice the rms residual against the published fraction is
  0.01, and the value does not decrease as the published fraction rises.
- The smooth ice fraction keeps a slope in liquid water, from a softer
  staggered weight still climbing at the published -0.8 edge and a
  Gaussian gate on a neighbour count of four. At equal wall time it
  is faster per CPU-hour than Q6 with a neighbour list and than
  ENVIRONMENTSIMILARITY to this hexagonal ice.
- The soft channel of the smooth ice fraction now stops at 0.65, so a
  liquid with no crystal can no longer outscore a seeded nucleus.
  Ideal ice Ic prints 0.99850519.

## [0.3.1] - 2026-09-02

### Changed

- Engine wrap pinned to seams-core v2.9.2.

## [0.3.0] - 2026-09-02

### Added

- SIGNATURE counts closed polyhedra by ring-size census on the same
  frame as the seeded ice score.
- GUESTS and GUEST_RADIUS place guest atoms at the periodic centroid
  of each SIGNATURE cage and report noccupied, nmultiple, and
  nfreeguest.
- LIBRARY accepts several key libraries at different hop counts; the
  deepest library that knows an atom names it.
- A PLUMED-NEST egg under nest/ with ice, brine, and hydrate inputs.
- Continuous integration builds the module and runs the driver tests
  on every push.

## [0.2.1] - 2026-09-02

### Changed

- Engine wrap pinned to the seams-core revision that ships with this
  cut.

## [0.2.0] - 2026-09-02

### Added

- HOPS (default 3) and LIBRARY on DSEAMS_CAGES.
- Components nclasses (distinct local keys on the mutual graph) and
  nnamed (molecules the library names).
- Driver test library_check.py: a cubic lattice, a library from
  seams fingerprint --emit-library, nnamed and nclasses on the
  driver.

### Changed

- Engine wrap follows seams-core with readcon-core v0.14.10.
- README.org documents install against an installed PLUMED, LOAD,
  fix plumed before fix npt, every keyword and component, the
  three driver tests, and the citation.

## [0.1.0] - 2026-09-02

### Added

- DSEAMS_CAGES, a PLUMED Colvar with no derivatives: cage occupancy
  (nice), largest connected cluster (nmax), cluster count (nclus),
  ice Ic / Ih / mixed (nic, nih, nmixed), CHILL+ (chillice, chillmax,
  chillinterfacial), and six-rings (sixrings).
- Keywords ATOMS, CUTOFF, CANDIDATE, K, LENGTH_SCALE, COMPLETE.
- IONS and ION_CUTOFF, with components nionice, nionfront, nionliq.
- Driver tests against walk_compare on the mW cubic fixture and the
  figshare nucleation deposit; ion lattice sites as ice.

### Fixed

- Component names omit underscores, which PLUMED reserves.
- Meson reads the PLUMED include directory from the last line of
  `plumed info --include-dir`.
