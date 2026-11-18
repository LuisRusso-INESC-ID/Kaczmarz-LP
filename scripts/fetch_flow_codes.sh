#!/bin/sh
# Download the reference max-flow codes used as extra baselines. They are not redistributed in
# this repository because of their licenses; `make` builds bin/flow-* from them if present.
#   HIPR 4.2     A. V. Goldberg, highest-label push-relabel (re-implementation of Cherkassky &
#                Goldberg's hi_pr). Academic, non-commercial license (LICENSE in the archive).
#                Source: Lehigh Max-Flow Challenge 2.0 reference solvers.
#   pseudo 3.23  B. Chandran and D. S. Hochbaum, pseudoflow (HPF). Academic license (see the
#                download page, riot.ieor.berkeley.edu/Applications/Pseudoflow/maxflow.html).
#   maxflow 3.04 Y. Boykov and V. Kolmogorov, BK augmenting paths with search trees. GPL.
# usage: scripts/fetch_flow_codes.sh [DIR]   (default $HOME/lpdeps/src/flow, the Makefile's FLOW)
set -eu
DIR=${1:-$HOME/lpdeps/src/flow}
mkdir -p "$DIR"
cd "$DIR"
get() {  # url file sha256
  [ -f "$2" ] || curl -fsSL -o "$2" "$1"
  echo "$3  $2" | sha256sum -c -
}
get https://coral.ise.lehigh.edu/wp-content/uploads/flows20/hipr42.tgz hipr42.tgz \
  13bd475c0c063f57e19652b264548030ac7e5c24125005080b857809b62866f4
get https://riot.ieor.berkeley.edu/Applications/Pseudoflow/pseudo-max-3.23_2.tar pseudo-max-3.23_2.tar \
  722bf7d8f8ad7dedb776171e71742e8f5643d1f8cb266b900e2147ff2b390903
get https://pub.ist.ac.at/~vnk/software/maxflow-v3.04.src.zip maxflow-v3.04.src.zip \
  1edf646f40674b2e89e113067ea825976eaff61c5295958820132895e4c6258b
rm -rf hipr42 pseudo_max maxflow-v3.04.src
tar -xzf hipr42.tgz
T=$(mktemp -d); tar -xf pseudo-max-3.23_2.tar -C "$T"
mv "$(find "$T" -type d -name pseudo_max | head -1)" pseudo_max; rm -rf "$T"
python3 -I -m zipfile -e maxflow-v3.04.src.zip .
ls -d "$DIR"/hipr42/main.c "$DIR"/pseudo_max/src/3.23/pseudo.c "$DIR"/maxflow-v3.04.src/maxflow.cpp
