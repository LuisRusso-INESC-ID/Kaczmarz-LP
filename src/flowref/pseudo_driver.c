/* pseudo_driver.c -- lpgraph driver for Chandran and Hochbaum's pseudoflow code, version 3.23
 * (academic license, fetched by scripts/fetch_flow_codes.sh and not distributed here).
 *
 * usage: flow-pseudoflow instance.max
 * Includes pseudo.c with its main() renamed, so the code is used unchanged; compiled as the
 * lowest-label FIFO variant (-DLOWEST_LABEL -DFIFO_BUCKET), which the authors report as
 * generally the fastest. t_solve covers initialization + phase 1, after which the minimum cut,
 * and so the max-flow value, is known (computed below as in pseudo.c's checkOptimality); t_flow
 * adds recoverFlow. Reading is excluded. Output: RESULT/HWM lines as bin/classic. */
#define main pseudo_original_main
#include "pseudo.c"
#undef main
#include <string.h>
#include <time.h>

static double now(void) {
  struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t);
  return t.tv_sec + 1e-9 * t.tv_nsec;
}
static long vm_hwm_kb(void) {
  char l[256]; long v = -1; FILE *f = fopen("/proc/self/status", "r");
  if (!f) return -1;
  while (fgets(l, sizeof l, f)) if (!strncmp(l, "VmHWM:", 6)) { v = atol(l + 6); break; }
  fclose(f); return v;
}

int main(int argc, char **argv) {
  double t0, t1, tcut, tflow;
  uint gap, i;
  ullint mincut = 0;
  if (argc < 2) { fprintf(stderr, "usage: flow-pseudoflow instance.max\n"); return 1; }
  if (!freopen(argv[1], "r", stdin)) { perror(argv[1]); return 1; }
  t0 = now();
  readDimacsFileCreateList();
  t1 = now();
  simpleInitialization();
  pseudoflowPhase1();
  tcut = now() - t1;
#ifdef LOWEST_LABEL
  gap = lowestStrongLabel;
#else
  gap = numNodes;
#endif
  for (i = 0; i < numArcs; ++i)
    if (arcList[i].from->label >= gap && arcList[i].to->label < gap) mincut += arcList[i].capacity;
  recoverFlow(gap);
  tflow = now() - t1;
  printf("RESULT status=ok obj=%llu t_solve=%.6f t_flow=%.6f t_total=%.6f\n", mincut, tcut, tflow, now() - t0);
  printf("HWM hwm_end_kb=%ld\n", vm_hwm_kb());
  return 0;
}
