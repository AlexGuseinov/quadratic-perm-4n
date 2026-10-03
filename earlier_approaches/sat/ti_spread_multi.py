"""Largest k such that some line spread of PG(n-1, 2) is totally isotropic for a
k-dimensional space W of alternating forms containing a form of rank `rank`.

B_1 is fixed to the canonical form of rank `rank`; B_2..B_k are variables; W must
have dimension k (every nonzero combination is a nonzero form). A chosen line
<u, v> must satisfy B_i(u, v) = 0 for all i.
Usage: python3 ti_spread_multi.py n rank k
"""

import os
import subprocess
import sys
from itertools import combinations

from ti_spread import B

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    n, rank, k = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    r = rank // 2
    N = 1 << n
    idx = list(combinations(range(n), 2))
    lines, seen = [], set()
    for x in range(1, N):
        for y in range(x + 1, N):
            if B(x, y, r) == 0:
                L = frozenset((x, y, x ^ y))
                if L not in seen:
                    seen.add(L)
                    lines.append((L, x, y))
    nv = len(lines)
    cls = []

    def new():
        nonlocal nv
        nv += 1
        return nv

    def xor(lits):
        acc = None
        for l in lits:
            if acc is None:
                acc = l
                continue
            o = new()
            cls.extend([[-o, acc, l], [-o, -acc, -l], [o, -acc, l], [o, acc, -l]])
            acc = o
        return acc  # None means constant 0

    # form variables for B_2..B_k
    bv = [[new() for _ in idx] for _ in range(k - 1)]
    canon = [1 if (i % 2 == 0 and j == i + 1 and i < rank) else 0 for (i, j) in idx]
    # exact cover
    through = {p: [] for p in range(1, N)}
    for t, (L, x, y) in enumerate(lines):
        for p in L:
            through[p].append(t + 1)
    for p in range(1, N):
        ls = through[p]
        cls.append(ls[:])
        m = len(ls)
        s = [new() for _ in range(m - 1)]
        cls.append([-ls[0], s[0]])
        for i in range(1, m - 1):
            cls += [[-ls[i], s[i]], [-s[i - 1], s[i]], [-ls[i], -s[i - 1]]]
        cls.append([-ls[m - 1], -s[m - 2]])
    # chosen line => t.i. for every variable form
    for t, (L, x, y) in enumerate(lines):
        coef = [((x >> i) & (y >> j) & 1) ^ ((x >> j) & (y >> i) & 1) for (i, j) in idx]
        for f in bv:
            val = xor([f[q] for q in range(len(idx)) if coef[q]])
            if val is not None:
                cls.append([-(t + 1), -val])
    # dim W = k: every nonzero combination (with or without B_1) is nonzero
    for c in range(1, 1 << (k - 1)):
        fs = [bv[i] for i in range(k - 1) if c >> i & 1]
        for with_b1 in (0, 1):
            coords = []
            for q in range(len(idx)):
                lit = xor([f[q] for f in fs])
                if with_b1 and canon[q]:
                    # coordinate = lit xor 1: nonzero iff lit false
                    coords.append(-lit if lit is not None else None)
                    if lit is None:
                        coords = ["TRUE"]
                        break
                else:
                    if lit is not None:
                        coords.append(lit)
            if "TRUE" in coords:
                continue
            cls.append([l for l in coords if l is not None])
    path = os.path.join(HERE, "cnf", f"ti_multi_n{n}_r{rank}_k{k}.cnf")
    with open(path, "w") as f:
        f.write(f"p cnf {nv} {len(cls)}\n")
        for c in cls:
            f.write(" ".join(map(str, c)) + " 0\n")
    out = subprocess.run([os.path.join(HERE, "..", "bin", "kissat"), "-q", path], capture_output=True, text=True).stdout
    st = [l for l in out.splitlines() if l.startswith("s ")]
    print(f"n={n} rank={rank} k={k}: vars={nv} clauses={len(cls)} -> {st[0] if st else 'no result'}", flush=True)


if __name__ == "__main__":
    main()
