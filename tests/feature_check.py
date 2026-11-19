#!/usr/bin/env python3
"""Query the structure-reader feature flags DSEAMS_CAGES exports.

  feature_check.py PLUMED_BIN MODULE_SO [--require-chemfiles] [--require-readcon]

DSEAMS_CAGES never opens a structure file itself (frames come from the MD
engine through PLUMED's getPosition()/getBox()), so there is nothing here to
read a .con, PDB or GRO file through: this test exists to prove the feature
query the module exports (dseams_plumed_has_chemfiles /
dseams_plumed_has_readcon, see src/DseamsCages.cpp) correctly reports what
the linked seams-core engine supports, matching the SEAMS_HAS_CHEMFILES /
SEAMS_HAS_READCON compile_args meson.build picked up from the engine
dependency. See README.org, Install.

libdseams_plumed.so has unresolved PLMD::* symbols (the PLUMED kernel
resolves them when LOAD dlopens the module; see meson.build), so a bare
ctypes.CDLL on it fails with "undefined symbol" the way any process that
has not already loaded libplumedKernel would. This script preloads the
kernel next to PLUMED_BIN with RTLD_GLOBAL first, the same resolution a
running `plumed driver` gives LOAD.

With --require-chemfiles / --require-readcon, the check fails (not skips)
when that reader is absent, so CI cannot silently pass without the engine
it means to test. Without those flags it only reports and never skips.
"""
from __future__ import annotations

import ctypes
import pathlib
import subprocess
import sys


def plumed_info(plumed: str, flag: str) -> str:
    out = subprocess.run([plumed, "info", flag], check=True, capture_output=True, text=True).stdout
    return out.strip().splitlines()[-1].strip()


def kernel_path(plumed: str) -> pathlib.Path:
    root = pathlib.Path(plumed_info(plumed, "--root"))
    soext = plumed_info(plumed, "--soext")
    # plumed info --root prints .../lib/plumed; the kernel sits beside it in
    # .../lib, the same place `plumed` itself prints as LibraryPath at
    # startup ("PLUMED: LibraryPath: .../lib/libplumedKernel.so").
    return root.parent / f"libplumedKernel.{soext}"


def main() -> int:
    args = sys.argv[1:]
    require_chemfiles = "--require-chemfiles" in args
    require_readcon = "--require-readcon" in args
    positional = [a for a in args if not a.startswith("--")]
    if len(positional) != 2:
        print(__doc__, file=sys.stderr)
        return 2
    plumed, module = positional

    kernel = kernel_path(plumed)
    if not kernel.is_file():
        print(f"no PLUMED kernel at {kernel}", file=sys.stderr)
        return 1
    ctypes.CDLL(str(kernel), mode=ctypes.RTLD_GLOBAL)

    lib = ctypes.CDLL(module)
    lib.dseams_plumed_has_chemfiles.restype = ctypes.c_bool
    lib.dseams_plumed_has_readcon.restype = ctypes.c_bool
    has_chemfiles = bool(lib.dseams_plumed_has_chemfiles())
    has_readcon = bool(lib.dseams_plumed_has_readcon())
    print(f"has_chemfiles={has_chemfiles} has_readcon={has_readcon}")

    if require_chemfiles and not has_chemfiles:
        print("chemfiles support was required but the linked engine lacks it", file=sys.stderr)
        return 1
    if require_readcon and not has_readcon:
        print("readcon support was required but the linked engine lacks it", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
