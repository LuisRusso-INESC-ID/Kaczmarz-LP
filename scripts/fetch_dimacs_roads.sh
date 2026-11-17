#!/bin/sh
# Download 9th DIMACS Challenge road networks (distance graphs) into data/.
# usage: scripts/fetch_dimacs_roads.sh [NY BAY COL FLA ...]   (default: NY BAY COL)
# Override the source with DIMACS_URL=... (directory holding USA-road-d.XXX.gr.gz).
# NOTE (2026-11-16): the challenge site returns 404; use scripts/fetch_zenodo_roads.sh (travel-time graphs).
set -e
URL=${DIMACS_URL:-https://www.diag.uniroma1.it/challenge9/data/USA-road-d}
cd "$(dirname "$0")/../data"
[ $# -gt 0 ] || set -- NY BAY COL
for r in "$@"; do
  f="USA-road-d.$r.gr"
  if [ ! -s "$f" ]; then
    curl -fsSL -o "$f.gz.part" "$URL/$f.gz" || { rm -f "$f.gz.part"; echo "download failed: $URL/$f.gz" >&2; exit 1; }
    gunzip -c "$f.gz.part" > "$f.part" && mv "$f.part" "$f" && rm -f "$f.gz.part"
  fi
  echo "$f: $(grep '^p ' "$f")"
done
# Run e.g.:  bin/classic dijkstra data/USA-road-d.NY.gr   /   bin/lpsolve highs-dual data/USA-road-d.NY.gr 600
