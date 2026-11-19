`libdseams_plumed.so` exports `dseams_plumed_has_chemfiles()` and
`dseams_plumed_has_readcon()`, reporting whether the linked seams-core
engine was built with its optional PDB/GRO/DCD (chemfiles) and `.con`
(readcon-core) readers. `DSEAMS_CAGES` logs the same two values at setup.
`DSEAMS_CAGES` never reads a structure file itself, so neither reader is
used by this action; see README.org, Install.
