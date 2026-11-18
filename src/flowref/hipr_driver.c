/* hipr_driver.c -- lpgraph driver for HIPR 4.2 (A. V. Goldberg; academic license, fetched by
 * scripts/fetch_flow_codes.sh and not distributed here).
 *
 * usage: flow-hipr instance.max
 * Includes HIPR's main.c with its main() renamed, so HIPR's code is used unchanged, and runs its
 * phases with a wall clock. t_solve covers init + stageOne, which yields the max-flow value
 * (a maximum preflow, as LEMON Preflow::runMinCut); t_flow adds stageTwo, which turns the
 * preflow into a flow. Reading is excluded. Output: RESULT/HWM lines as bin/classic. */
#define main hipr_original_main
#include "main.c"
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

int main(int argc, char *argv[]) {
  double t0, t1, tcut, tflow, value;
  if (argc < 2) { fprintf(stderr, "usage: flow-hipr instance.max\n"); return 1; }
  if (!freopen(argv[1], "r", stdin)) { perror(argv[1]); return 1; }
  t0 = now();
  globUpdtFreq = GLOB_UPDT_FREQ; globUpdtFreqOrig = globUpdtFreq;
  parse(&n, &m, &mInp, &nodes, &arcs, &cap, &source, &sink, &nMin);
  if (allocDS()) { printf("RESULT status=fail(alloc)\n"); return 1; }
  t1 = now();
  if (init()) stageOne(); else flow = 0;
  tcut = now() - t1;
  value = flow;
  stageTwo();
  tflow = now() - t1;
  printf("RESULT status=ok obj=%.0f t_solve=%.6f t_flow=%.6f t_total=%.6f\n", value, tcut, tflow, now() - t0);
  printf("HWM hwm_end_kb=%ld\n", vm_hwm_kb());
  return 0;
}
