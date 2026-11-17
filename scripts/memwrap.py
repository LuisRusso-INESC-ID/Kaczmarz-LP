#!/usr/bin/env python3
"""memwrap.py MEMLIMIT TIMEOUT cmd... -- run cmd, forward its stdout/stderr, and report
wall time, peak RSS (ru_maxrss of the child via wait4) and whether it was killed, as a
'MEMWRAP {json}' line on stderr.  LPBENCH_CPU=k pins the solver process to CPU k."""
import json, os, resource, signal, sys, time

limit, timeout, cmd = int(sys.argv[1]), float(sys.argv[2]), sys.argv[3:]
t0 = time.monotonic()
pid = os.fork()
if pid == 0:
    resource.setrlimit(resource.RLIMIT_AS, (limit, limit))
    if os.environ.get("LPBENCH_CPU"):   # pin the solver (not this wrapper) to one core
        os.sched_setaffinity(0, {int(os.environ["LPBENCH_CPU"])})
    os.execv(cmd[0], cmd)
killed = False
while True:
    wpid, status, ru = os.wait4(pid, os.WNOHANG)
    if wpid == pid:
        break
    if time.monotonic() - t0 > timeout:
        os.kill(pid, signal.SIGKILL); killed = True
        wpid, status, ru = os.wait4(pid, 0)
        break
    time.sleep(0.005)
wall = time.monotonic() - t0
rc = os.waitstatus_to_exitcode(status)
sys.stderr.write("MEMWRAP " + json.dumps({"wall": wall, "maxrss_kb": ru.ru_maxrss, "killed": killed, "returncode": rc}) + "\n")
