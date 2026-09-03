# Changelog

All notable changes to this project are documented in this file.

The format is [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Unreleased notes live in [`changelog.d/`](changelog.d/) and are assembled
by [towncrier](https://towncrier.readthedocs.io/).

<!-- towncrier release notes start -->

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
