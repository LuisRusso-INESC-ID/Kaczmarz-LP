#!/bin/sh
# Download 9th DIMACS Challenge road networks (distance graphs) into data/.
# Not reachable from the sandbox where the pilot ran; run this on your own machine.
# usage: scripts/fetch_dimacs_roads.sh NY BAY COL FLA
set -e
cd "$(dirname "$0")/../data"
for r in "${@:-NY BAY COL}"; do
  f="USA-road-d.$r.gr"
  [ -f "$f" ] || { curl -fsSL "https://www.diag.uniroma1.it/challenge9/data/USA-road-d/$f.gz" | gunzip > "$f"; }
  echo "$f: $(grep '^p ' "$f")"
done
# Run e.g.:  bin/classic dijkstra data/USA-road-d.NY.gr   /   bin/lpsolve highs-dual data/USA-road-d.NY.gr 600
