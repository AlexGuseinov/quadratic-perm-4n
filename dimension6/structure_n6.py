"""Structure of constant-rank-(n-2), dimension-n subspaces W of Alt(n, F_2)
with s_a = 2 for every a (the DU-4 case), explored on many n = 6 solutions.

For each solution we check:
  (1) W-isotropic lines {<u,v> : B(u,v) = 0 for all B in W} form a line spread;
  (2) the radicals R_b (b != 0) are exactly those spread lines;
  (3) each radical is shared by exactly 3 nonzero b forming {b1, b2, b1+b2};
  (4) the lines S_a = {b : a in R_b} form a line spread of W.
Solutions are sampled with blocking clauses on the basis variables.
"""

import random
import sys
from itertools import combinations

from pysat.solvers import Solver

import os, sys  # repository layout: shared modules are in ../core
sys.path.insert(1, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "core"))
from constrank_sat import build, decode, combo, rank_rows


def radical(M, n):
    return frozenset(a for a in range(1 << n)
                     if all(bin(M[i] & a).count("1") % 2 == 0 for i in range(n)))


def analyse(mats, n):
    d = len(mats)
    forms = {c: combo(mats, c, n) for c in range(1, 1 << d)}
    rads = {c: radical(M, n) for c, M in forms.items()}

    def B(M, u, v):
        acc = 0
        for i in range(n):
            if (u >> i) & 1:
                acc ^= bin(M[i] & v).count("1") & 1
        return acc

    iso = set()
    for u in range(1, 1 << n):
        for v in range(u + 1, 1 << n):
            if all(B(mats[k], u, v) == 0 for k in range(d)):
                iso.add(frozenset({0, u, v, u ^ v}))
    pts = [p for L in iso for p in L if p]
    spread = len(pts) == (1 << n) - 1 and len(set(pts)) == len(pts)
    rad_set = set(rads.values())
    radicals_are_spread = rad_set == iso
    by_rad = {}
    for c, R in rads.items():
        by_rad.setdefault(R, []).append(c)
    triples = all(len(v) == 3 and v[0] ^ v[1] == v[2] for v in
                  (sorted(x) for x in by_rad.values()))
    S = {}
    for a in range(1, 1 << n):
        S[a] = frozenset([0] + [c for c, R in rads.items() if a in R])
    s_lines = set(S.values())
    s_pts = [p for L in s_lines for p in L if p]
    s_spread = len(s_pts) == (1 << d) - 1 and len(set(s_pts)) == len(s_pts)
    return spread, radicals_are_spread, triples, s_spread


def main():
    n = 6
    samples = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    enc, basis, pairs = build(n, n - 2, n, du4=True)
    basis_vars = [v for row in basis for v in row]
    stats = {}
    with Solver(name="cadical195", bootstrap_with=enc.clauses) as s:
        for t in range(samples):
            random.seed(t)
            assumptions = [v if random.random() < 0.5 else -v
                           for v in random.sample(basis_vars[len(pairs):], 4)]
            if not s.solve(assumptions=assumptions):
                if not s.solve():
                    break
            model = s.get_model()
            mats = decode(model, basis, pairs, n)
            key = analyse(mats, n)
            stats[key] = stats.get(key, 0) + 1
            pos = set(x for x in model if x > 0)
            s.add_clause([-v if v in pos else v for v in basis_vars])
    print("(isotropic lines form spread, radicals = those lines, "
          "radicals shared in triples, S_a lines form spread): count")
    for k, v in stats.items():
        print(" ", k, v)


if __name__ == "__main__":
    main()
