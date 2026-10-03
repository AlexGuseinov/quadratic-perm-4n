"""Does F_2^n admit a line spread all of whose lines are totally isotropic for a
fixed alternating form B of rank 2r?  (B = e0^e1 + e2^e3 + ... , r terms; all
forms of the same rank are equivalent.)

Relevance: for a quadratic DU-4 permutation, the isotropic spread Sigma (L3) is
totally isotropic for every component form B_b, and every B_b has rank n-2 (L2).
So if no B-t.i. line spread exists for rank n-2, D(n) has no solution.

Exact cover of the 2^n - 1 points by B-t.i. lines; solved with kissat.
Usage: python3 ti_spread.py n rank [--kissat PATH]
"""

import os
import subprocess
import sys
from itertools import combinations


def B(x, y, r):
    s = 0
    for i in range(r):
        a0, a1 = (x >> 2 * i) & 1, (x >> 2 * i + 1) & 1
        b0, b1 = (y >> 2 * i) & 1, (y >> 2 * i + 1) & 1
        s ^= (a0 & b1) ^ (a1 & b0)
    return s

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    n, rank = int(sys.argv[1]), int(sys.argv[2])
    r = rank // 2
    kissat = sys.argv[sys.argv.index("--kissat") + 1] if "--kissat" in sys.argv else os.path.join(HERE, "..", "bin", "kissat")
    N = 1 << n
    lines = []
    seen = set()
    for x in range(1, N):
        for y in range(x + 1, N):
            if B(x, y, r) == 0:
                L = frozenset((x, y, x ^ y))
                if L not in seen:
                    seen.add(L)
                    lines.append(L)
    through = {p: [] for p in range(1, N)}
    for i, L in enumerate(lines):
        for p in L:
            through[p].append(i + 1)
    clauses = []
    nv = len(lines)
    for p in range(1, N):
        ls = through[p]
        clauses.append(ls[:])
        # at most one: sequential counter (Sinz)
        k = len(ls)
        if k > 1:
            s = list(range(nv + 1, nv + k))
            nv += k - 1
            clauses.append([-ls[0], s[0]])
            for i in range(1, k - 1):
                clauses.append([-ls[i], s[i]])
                clauses.append([-s[i - 1], s[i]])
                clauses.append([-ls[i], -s[i - 1]])
            clauses.append([-ls[k - 1], -s[k - 2]])
    # symmetry breaking is omitted (answers are small instances)
    os.makedirs(os.path.join(HERE, "cnf"), exist_ok=True)
    path = os.path.join(HERE, "cnf", f"ti_spread_n{n}_r{rank}.cnf")
    with open(path, "w") as f:
        f.write(f"p cnf {nv} {len(clauses)}\n")
        for cl in clauses:
            f.write(" ".join(map(str, cl)) + " 0\n")
    print(f"n={n} rank={rank}: t.i. lines={len(lines)} vars={nv} clauses={len(clauses)}", flush=True)
    out = subprocess.run([kissat, "-q", path], capture_output=True, text=True).stdout
    status = [l for l in out.splitlines() if l.startswith("s ")]
    print(" ", status[0] if status else out[:200])
    if status and "UNSAT" not in status[0]:
        vals = set()
        for l in out.splitlines():
            if l.startswith("v "):
                vals.update(int(t) for t in l[2:].split() if int(t) > 0)
        chosen = [lines[i - 1] for i in range(1, len(lines) + 1) if i in vals]
        pts = [p for L in chosen for p in L]
        ok = len(pts) == N - 1 and len(set(pts)) == N - 1 and all(
            B(a, b, r) == 0 for L in chosen for a in L for b in L)
        print(f"  spread with {len(chosen)} lines; independent check: {ok}")


if __name__ == "__main__":
    main()
