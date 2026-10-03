"""Write ONE cube of the D(8) SB3 model (base clauses + unit clauses) to a file.
Cube m of depth d fixes the first d free c-variables to the bits of m.
Usage (from earlier_approaches/sat): python3 vps/cube_cnf.py depth m out.cnf"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "core"))
from du4perm_sat import build
depth, m, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
enc, c, ell = build(8, sb=True, sb3=True)
free = [v for p in sorted(c) for v in c[p] if isinstance(v, int) and not isinstance(v, bool)]
units = [free[i] if (m >> i) & 1 else -free[i] for i in range(depth)]
with open(out, "w") as f:
    f.write(f"c D(8) sb3 cube depth={depth} m={m} units={units}\n")
    f.write(f"p cnf {enc.pool.top} {len(enc.clauses) + depth}\n")
    for cl in enc.clauses:
        f.write(" ".join(map(str, cl)) + " 0\n")
    for u in units:
        f.write(f"{u} 0\n")
