#!/bin/bash
# Build the LP solvers and LEMON used by lpgraph from source into $PREFIX (default ./deps/install).
# Mirrors what was done for the pilot (HiGHS master ~1.15, Clp/CoinUtils master, GLPK, LEMON 1.3.1
# headers, OR-Tools C++ with GLOP and PDLP).  Then build lpgraph with
#     make DEPS=$PREFIX LEMON=$DEPSRC/lemon/include
# Usage: scripts/build_deps.sh [PREFIX] [SRCDIR]
set -euo pipefail
PREFIX=${1:-$PWD/deps/install}
SRC=${2:-$PWD/deps}
JOBS=${JOBS:-$(nproc)}
mkdir -p "$SRC" "$PREFIX"
cd "$SRC"
clone() { [ -d "$(basename "$1")" ] || git clone --depth 1 "https://github.com/$1.git"; }

clone ERGO-Code/HiGHS
clone coin-or/CoinUtils
clone coin-or/Clp
clone firedrakeproject/glpk      # GitHub mirror of GNU GLPK; the GNU tarball works as well
clone seqan/lemon                # LEMON 1.3.1 headers (header-only use)
clone google/or-tools

# HiGHS
cmake -S HiGHS -B HiGHS/build -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX="$PREFIX"
cmake --build HiGHS/build -j"$JOBS" && cmake --install HiGHS/build
# GLPK
(cd glpk && ./configure --prefix="$PREFIX" --disable-static CFLAGS=-O2 && make -j"$JOBS" && make install)
# CoinUtils + Clp (no BLAS/LAPACK, as in the pilot)
(cd CoinUtils && ./configure --prefix="$PREFIX" --disable-static --without-lapack --without-blas CXXFLAGS=-O2 && make -j"$JOBS" && make install)
(cd Clp && PKG_CONFIG_PATH="$PREFIX/lib/pkgconfig" ./configure --prefix="$PREFIX" --disable-static --without-lapack --without-blas CXXFLAGS=-O2 && make -j"$JOBS" && make install)
# OR-Tools (C++ only; GLOP and PDLP are always built).  BUILD_DEPS fetches abseil, protobuf,
# re2 and Eigen.  In the pilot sandbox gitlab.com was unreachable, so Eigen was fetched from a
# GitHub mirror and system zlib/bzip2 were used; on a normal machine the defaults below work.
cmake -S or-tools -B or-tools/build -DCMAKE_BUILD_TYPE=Release -DBUILD_DEPS=ON -DBUILD_CXX=ON \
  -DBUILD_PYTHON=OFF -DBUILD_SAMPLES=OFF -DBUILD_EXAMPLES=OFF -DBUILD_TESTING=OFF \
  -DUSE_SCIP=OFF -DUSE_COINOR=OFF -DUSE_HIGHS=OFF -DUSE_GLPK=OFF \
  -DCMAKE_INSTALL_PREFIX="$PREFIX/ortools"
cmake --build or-tools/build -j"$JOBS" && cmake --install or-tools/build
echo "done: make DEPS=$PREFIX LEMON=$SRC/lemon/include"
