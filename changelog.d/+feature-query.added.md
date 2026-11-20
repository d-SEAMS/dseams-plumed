`libdseams_plumed.so` exports `dseams_plumed_has_chemfiles()` and
`dseams_plumed_has_readcon()`. `tests/feature_check.py` checks those two
flags. When the seams-core wrap is built with both readers,
`tests/reader_fail_check.py` runs `seams read` on a nonexistent file, an
empty file, and a mode-000 file for each of PDB, GRO, DCD, and `.con`,
and requires exit 2, empty stdout, and a stderr line `cannot read <path>`.
`DSEAMS_CAGES` does not read those files.
