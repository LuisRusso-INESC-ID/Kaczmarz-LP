#!/bin/sh
# Hop-depth experiment (results/depth/): grids with 4,096 nodes and lengths U[1,1000] but aspect
# ratios 64x64 ... 4x1024, so the shortest-path-tree depth grows ~10x at fixed size; own PDHG and
# OR-Tools PDLP to tolerance 1e-6 (900 s limit), pinned to CPU ${LPBENCH_CPU:-2}; Dijkstra reference.
# Also the hop depths of the benchmark instances that the first-order methods reach.
set -e
ROOT=$(cd "$(dirname "$0")/.." && pwd); cd "$ROOT"
CPU=${LPBENCH_CPU:-2}; OUT=results/depth; mkdir -p $OUT data
: > $OUT/runs.txt; : > $OUT/ref.txt
for g in 64x64 32x128 16x256 8x512 4x1024; do
  f=data/depth-grid-$g.gr; bin/gen sp-grid ${g%x*} ${g#*x} 1 > $f
  echo "$g $(bin/classic dijkstra $f 2>/dev/null | tr ' ' '\n' | grep '^obj=')" >> $OUT/ref.txt
  for m in own-pdhg ortools-pdlp; do
    b=bin/lpsolve; [ $m = ortools-pdlp ] && b=bin/lpsolve-ort
    echo "$g $m $(taskset -c $CPU $b $m $f 900 1e-6 2>/dev/null | tr ' ' '\n' | grep -E '^(status|iters|t_solve|obj)=' | tr '\n' ' ')" >> $OUT/runs.txt
  done
done
python3 -I scripts/hop_depth.py data/depth-grid-*.gr > $OUT/depths.txt
python3 -I scripts/hop_depth.py data/sp-grid-32x32-s1.gr data/sp-grid-64x64-s1.gr data/sp-grid-128x128-s1.gr \
  data/sp-road-32x32-s1.gr data/sp-road-64x64-s1.gr data/sp-road-128x128-s1.gr data/sp-rand-4000x4-s1.gr \
  data/sp-rand-16000x4-s1.gr data/sp-rand-64000x4-s1.gr data/USA-road-t.NW-40k.gr data/USA-road-t.NW-160k.gr >> $OUT/depths.txt
