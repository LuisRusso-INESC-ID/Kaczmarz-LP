#!/bin/bash
# Build the LP solvers and LEMON used by lpgraph from source into $PREFIX (default ./deps/install).
# Mirrors what was used for the paper's runs (HiGHS 1.16, Clp/CoinUtils master with SuiteSparse
# AMD+CHOLMOD for clp-barrier-chol, GLPK 5.0, LEMON 1.3.1 headers, OR-Tools 9.15 C++ with GLOP and
# PDLP).  The reference max-flow codes come separately: scripts/fetch_flow_codes.sh.  Then build with
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
clone seqan/lemon                # LEMON 1.3.1 headers (header-only use)
clone google/or-tools

# HiGHS
cmake -S HiGHS -B HiGHS/build -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX="$PREFIX"
cmake --build HiGHS/build -j"$JOBS" && cmake --install HiGHS/build
# GLPK 5.0 from GNU (-std=gnu17: GCC 15+ defaults to C23, where GLPK's minisat typedef of bool fails)
[ -f glpk-5.0.tar.gz ] || curl -fsSLO https://ftp.gnu.org/gnu/glpk/glpk-5.0.tar.gz
echo "4a1013eebb50f728fc601bdd833b0b2870333c3b3e5a816eeba921d95bec6f15  glpk-5.0.tar.gz" | sha256sum -c -
tar -xzf glpk-5.0.tar.gz
(cd glpk-5.0 && ./configure --prefix="$PREFIX" --disable-static CFLAGS="-O2 -std=gnu17" && make -j"$JOBS" && make install)
# SuiteSparse 7.14.1: AMD and CHOLMOD only, on the system's reference BLAS/LAPACK, no OpenMP
[ -f SuiteSparse-7.14.1.tar.gz ] || curl -fsSL -o SuiteSparse-7.14.1.tar.gz \
  https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/refs/tags/v7.14.1.tar.gz
echo "81e560e1f74546df139edb765b3f5bc865866da23062312ffe8fd821063c8397  SuiteSparse-7.14.1.tar.gz" | sha256sum -c -
tar -xzf SuiteSparse-7.14.1.tar.gz
cmake -S SuiteSparse-7.14.1 -B SuiteSparse-7.14.1/build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DSUITESPARSE_ENABLE_PROJECTS="suitesparse_config;amd;camd;colamd;ccolamd;cholmod" -DSUITESPARSE_USE_CUDA=OFF \
  -DSUITESPARSE_USE_OPENMP=OFF -DCHOLMOD_USE_OPENMP=OFF -DBUILD_STATIC_LIBS=OFF -DSUITESPARSE_DEMOS=OFF -DBUILD_TESTING=OFF
cmake --build SuiteSparse-7.14.1/build -j"$JOBS" && cmake --install SuiteSparse-7.14.1/build
# CoinUtils + Clp (no BLAS/LAPACK for Clp itself; AMD+CHOLMOD for the clp-barrier-chol backend;
# Clp's default barrier, clp-barrier, still uses its native Cholesky)
(cd CoinUtils && ./configure --prefix="$PREFIX" --disable-static --without-lapack --without-blas CXXFLAGS=-O2 && make -j"$JOBS" && make install)
SS="-L$PREFIX/lib -Wl,-rpath,$PREFIX/lib"
(cd Clp && PKG_CONFIG_PATH="$PREFIX/lib/pkgconfig" ./configure --prefix="$PREFIX" --disable-static --without-lapack --without-blas CXXFLAGS=-O2 \
   --with-cholmod-lflags="$SS -lcholmod -lamd -lcolamd -lsuitesparseconfig" --with-cholmod-cflags="-I$PREFIX/include/suitesparse" \
   --with-amd-lflags="$SS -lamd -lsuitesparseconfig" --with-amd-cflags="-I$PREFIX/include/suitesparse" && make -j"$JOBS" && make install)
# OR-Tools (C++ only; GLOP and PDLP are always built).  BUILD_DEPS fetches abseil, protobuf,
# re2 and Eigen.  In the pilot sandbox gitlab.com was unreachable, so Eigen was fetched from a
# GitHub mirror and system zlib/bzip2 were used; on a normal machine the defaults below work.
cmake -S or-tools -B or-tools/build -DCMAKE_BUILD_TYPE=Release -DBUILD_DEPS=ON -DBUILD_CXX=ON \
  -DBUILD_PYTHON=OFF -DBUILD_SAMPLES=OFF -DBUILD_EXAMPLES=OFF -DBUILD_TESTING=OFF \
  -DUSE_SCIP=OFF -DUSE_COINOR=OFF -DUSE_HIGHS=OFF -DUSE_GLPK=OFF \
  -DCMAKE_INSTALL_PREFIX="$PREFIX/ortools"
cmake --build or-tools/build -j"$JOBS" && cmake --install or-tools/build
echo "done: make DEPS=$PREFIX LEMON=$SRC/lemon/include"
