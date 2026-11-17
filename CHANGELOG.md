# Changelog

All notable changes to this project are documented in this file.

The format is [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Unreleased notes live in [`changelog.d/`](changelog.d/) and are assembled
by [towncrier](https://towncrier.readthedocs.io/).

<!-- towncrier release notes start -->

## [0.4.1] - 2026-11-17

### Fixed

- The PLUMED-NEST archive now contains everything needed to run. Its source compiles with `plumed mklib` under gcc and libstdc++ and under clang and libc++. It keeps the upstream identifiers as they are. The archive bundles only the Eigen headers the source includes, together with the licences they require. It needs PLUMED 2.10 or later.


## [0.4.0] - 2026-11-17

### Added

- SIGNATURE accepts the named tables seams-core added after the v2.9.2
  engine in v0.3.1: 51264, 51268, sh, and sH, beside sodalite, alpha,
  512, 51262, hc, and ddc.

### Changed

- Engine wrap pins the seams-core tag v2.14.0. v0.3.1 wrapped engine
  v2.9.2. The tag includes the unreleased pins since then: v2.10.0 (F4,
  Zeron, TUM offload, and the water-type mask), ice XXI fingerprints,
  the audit fixes for undefined behavior and the IRA residual, the j815
  leftover batch, CHILL+ device-safe target loops, and the RDF
  minimum-image correction. The module reads the library field colored.

### Fixed

- CHILL+ counts follow the seams-core correction on a cutoff list. The
  old count used bonds after the neighbour row had already been cut to
  four. A crystal row needs four neighbours inside the cutoff.
- LOAD examples pass `./libdseams_plumed.so`. The dynamic loader does
  not search the working directory for a bare filename.


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
