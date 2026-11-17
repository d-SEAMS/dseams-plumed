#!/usr/bin/env bash
# Build PLUMED from source with clang and libc++ so that LOAD FILE compiles
# the archive with the same standard library as the kernel it loads into.
# Usage: build_plumed_libcxx.sh ENV_PREFIX PLUMED_SOURCE INSTALL_PREFIX
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "usage: $0 ENV_PREFIX PLUMED_SOURCE INSTALL_PREFIX" >&2
  exit 2
fi
env_prefix="$(realpath "$1")"
source_dir="$(realpath "$2")"
mkdir -p "$3"
install_prefix="$(realpath "$3")"

export PATH="$env_prefix/bin:$PATH"
cd "$source_dir"
./configure --prefix="$install_prefix" \
  CC=clang CXX="clang++ -stdlib=libc++" \
  --enable-rpath --disable-mpi --disable-openmp --disable-python \
  --disable-doc --disable-libsearch
make -j"$(nproc)"
make install
