"""Run the exact D(8) model on one half of the search space.
Half h in {0,1}: assumption on the first entry of c_12 (a free variable)."""
import sys, time
from pysat.solvers import Solver
import os, sys  # repository layout: shared modules are in ../../core
sys.path.insert(1, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "core"))
from du4perm_sat import build, extract, lut, check
h = int(sys.argv[1])
n = 8
t0 = time.time()
enc, c, ell = build(n)
v = c[(1, 2)][0]
print(f"D(8) half {h}: vars={enc.pool.top} clauses={len(enc.clauses)} "
      f"assume {'+' if h else '-'}{v} (build {time.time()-t0:.1f}s)", flush=True)
with Solver(name="cadical195", bootstrap_with=enc.clauses) as s:
    t1 = time.time()
    res = s.solve(assumptions=[v if h else -v])
    dt = time.time() - t1
    if res:
        C, L = extract(s.get_model(), c, ell, n)
        table = lut(C, L, n)
        perm, du = check(table, n)
        print(f"SAT in {dt:.1f}s: permutation={perm} DU={du}", flush=True)
        print("LUT:", table, flush=True)
        open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "logs", f"SOLUTION_D8_half{h}.txt"), "w").write(str(table))
    else:
        print(f"UNSAT in {dt:.1f}s", flush=True)
