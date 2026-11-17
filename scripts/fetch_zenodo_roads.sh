#!/bin/sh
# Download the 9th DIMACS Challenge USA road graphs (travel time, largest connected component) from
# Zenodo record 10.5281/zenodo.18026780 (T. Abeywickrama, 2025; 635 MB) and convert them to DIMACS
# .gr/.co in data/:
#   USA-road-t.NW.gr        1,089,933 nodes,  2,545,844 arcs
#   USA-road-t.NW-{40k,160k,640k}.gr   BFS regions of NW around its centre (scripts/convert_roads.py)
#   USA-road-t.USA.gr      23,947,347 nodes, 57,708,624 arcs   (only with ALL=1)
# The original challenge site (diag.uniroma1.it/challenge9) returns 404 since 2026.
# usage: [ALL=1] [PY=python3-with-numpy] scripts/fetch_zenodo_roads.sh
set -e
ROOT=$(cd "$(dirname "$0")/.." && pwd)
PY=${PY:-python3}
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
mkdir -p "$ROOT/data"
curl -fL -o "$TMP/usa.tar.gz" https://zenodo.org/api/records/18026780/files/USA_road_network_and_poi.tar.gz/content
tar -xzf "$TMP/usa.tar.gz" -C "$TMP" data/NW-t.graph data/NW-t.coordinates
"$PY" -I "$ROOT/scripts/convert_roads.py" convert "$TMP/data/NW-t.graph" "$TMP/data/NW-t.coordinates" "$ROOT/data/USA-road-t.NW"
"$PY" -I "$ROOT/scripts/convert_roads.py" cut "$ROOT/data/USA-road-t.NW" 40000 160000 640000
if [ "${ALL:-0}" = 1 ]; then
  tar -xzf "$TMP/usa.tar.gz" -C "$TMP" data/USA-t.graph data/USA-t.coordinates
  # streaming conversion (the graph is 1.3 GB of text)
  awk 'NR==1{print "c 9th DIMACS Challenge USA road graph USA (travel time),"; print "c largest connected component, via Zenodo 10.5281/zenodo.18026780"; print "p sp",$3,$4; next} {print "a",$1+1,$2+1,$3}' \
    "$TMP/data/USA-t.graph" > "$ROOT/data/USA-road-t.USA.gr"
  awk 'NR==1{print "p aux sp co",$3; next}{print "v",$1+1,$2,$3}' "$TMP/data/USA-t.coordinates" > "$ROOT/data/USA-road-t.USA.co"
fi
ls -la "$ROOT"/data/USA-road-t.*
