"""Write the exact D(n) SAT model (du4perm_sat.build) as DIMACS.
Optional unit assumptions: --assume LIT[,LIT...] (added as unit clauses).
Usage: python3 export_dimacs.py n out.cnf [--assume 57] [--nosb]"""
import sys
import os, sys  # repository layout: shared modules are in ../../core
sys.path.insert(1, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "core"))
from du4perm_sat import build
n, out = int(sys.argv[1]), sys.argv[2]
units = []
if "--assume" in sys.argv:
    units = [int(x) for x in sys.argv[sys.argv.index("--assume") + 1].split(",")]
enc, c, ell = build(n, sb="--nosb" not in sys.argv, sb3="--sb3" in sys.argv, perm="--noperm" not in sys.argv, desarg="--desarg" in sys.argv, closed="--closed" in sys.argv)
cls = enc.clauses + [[u] for u in units]
with open(out, "w") as f:
    f.write(f"c D({n}) quadratic DU-4 permutation model, sb={'--nosb' not in sys.argv}, sb3={'--sb3' in sys.argv}, perm={'--noperm' not in sys.argv}, desarg={'--desarg' in sys.argv}, closed={'--closed' in sys.argv}, units={units}\n")
    f.write(f"p cnf {enc.pool.top} {len(cls)}\n")
    for cl in cls:
        f.write(" ".join(map(str, cl)) + " 0\n")
print(out, enc.pool.top, len(cls))
