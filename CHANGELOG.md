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

### Changed

- DSEAMS_SMOOTH_CHILL evaluates on PLUMED's NUM_THREADS. The pair list
  order is unchanged, so one thread and several threads return the same
  value, forces and virial.
- Eight walkers and 28 ns of well-tempered metadynamics at 210 K put
  the slab free-energy minimum at smooth fraction 0.930. At equal wall
  time the smooth fraction reaches a published ice fraction of 0.8 and
  Q6 MEAN does not.
- Engine wrap follows HaoZeke/seams-core 988d3996 (published CHILL+ rule
  and the ovito-chill+ alternative).
- Engine wrap follows seams-core 67628524 (v2.10.0 plus ice XXI).
- Engine wrap follows seams-core 91235ec1 (CHILL+ device-safe target loops).
- Engine wrap follows seams-core bdc9b74 (RDF MIC-once and ice XXI fingerprint).
- Engine wrap follows seams-core eb1a4985 (j815 leftover batch).
- Engine wrap follows seams-core fa1c8a73 (audit UB and IRA residual).
- Engine wrap follows the seams-core 2.14.0 release (628b0af9).
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
- The smooth ice fraction keeps a slope in liquid water. A softer
  staggered weight, still climbing at the published -0.8 edge, and a
  Gaussian gate on a neighbour count of four leave ideal ice Ic at
  0.99544. Eight walkers at 210 K, one started from the liquid, put
  the free-energy minimum at fraction 0.935. At equal wall time it
  is faster per CPU-hour than Q6 with a neighbour list and than
  ENVIRONMENTSIMILARITY to this hexagonal ice.
- The soft channel of the smooth ice fraction now stops at 0.65, so a
  liquid with no crystal can no longer outscore a seeded nucleus.
  Ideal ice Ic prints 0.99850519.


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
