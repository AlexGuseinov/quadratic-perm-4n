"""Reproduces the n = 6 observations in NOTES_lemmas.md ("Observations"):

  (a) Pluecker span of the isotropic spread = C(6,2) - 6 = 9, i.e. Sigma^perp = W;
  (b) self-adjoint algebra of W has dimension 1 (only scalars); x^5 gives 2 (= F_4);
  (c) number of "closed" pairs of spread lines: 50 for the samples, 210 for x^5;
  (d) Desarguesian spread for n = 8: dim Sigma^perp = 12.

Samples: random quadratic DU-4 permutations of F_2^6 from the exact SAT model.
Usage: python3 observations_n6.py [samples]
"""

import random
import sys
from collections import Counter
from itertools import combinations

from pysat.solvers import Solver

import os, sys  # repository layout: shared modules are in ../core
sys.path.insert(1, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "core"))
from du4perm_sat import build, extract, lut
from quadtools import power_map, rank_rows


def spread_of(tab, n):
    N = 1 << n
    beta = lambda u, v: tab[u ^ v] ^ tab[u] ^ tab[v] ^ tab[0]
    return list({frozenset(v for v in range(N) if beta(a, v) == 0) for a in range(1, N)})


def pluecker(u, v, idx):
    m = 0
    for k, (i, j) in enumerate(idx):
        if ((u >> i) & (v >> j) & 1) ^ ((u >> j) & (v >> i) & 1):
            m |= 1 << k
    return m


def span_dim(S, n):
    idx = list(combinations(range(n), 2))
    vecs = [pluecker(*sorted(x for x in L if x)[:2], idx) for L in S]
    return rank_rows(vecs), len(idx)


def selfadjoint_dim(tab, n):
    beta = lambda u, v: tab[u ^ v] ^ tab[u] ^ tab[v] ^ tab[0]
    rows = []
    for k in range(n):
        M = [[(beta(1 << i, 1 << j) >> k) & 1 for j in range(n)] for i in range(n)]
        for i in range(n):
            for j in range(n):
                r = 0
                for t in range(n):
                    if M[t][j]:
                        r ^= 1 << (t * n + i)
                    if M[i][t]:
                        r ^= 1 << (t * n + j)
                rows.append(r)
    return n * n - rank_rows(rows)


def closed_pairs(S):
    pt = {p: i for i, L in enumerate(S) for p in L if p}
    c = 0
    for A, B in combinations(S, 2):
        U = {a ^ b for a in A for b in B}
        if all(S[pt[p]] <= U for p in U if p):
            c += 1
    return c


def omega(x, n):
    y = 0
    for t in range(0, n, 2):
        a, b = (x >> t) & 1, (x >> (t + 1)) & 1
        y |= (b << t) | ((a ^ b) << (t + 1))
    return y


def main():
    samples = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    n = 6
    x5 = power_map(5, 6, 0b1000011)
    S = spread_of(x5, n)
    d, tot = span_dim(S, n)
    print(f"x^5 on GF(64): span={d}, dim Sigma^perp={tot - d}, "
          f"self-adjoint dim={selfadjoint_dim(x5, n)}, closed pairs={closed_pairs(S)}")
    enc, c, ell = build(n, sb=True, sb3=True)
    free = [v for p, vs in c.items() for v in vs if isinstance(v, int) and not isinstance(v, bool)]
    st = Counter()
    with Solver(name="cadical195", bootstrap_with=enc.clauses) as s:
        for t in range(samples):
            random.seed(20000 + t)
            assum = [v if random.random() < 0.5 else -v for v in random.sample(free, 8)]
            if not s.solve(assumptions=assum) and not s.solve():
                break
            m = s.get_model()
            C, L = extract(m, c, ell, n)
            tab = lut(C, L, n)
            S = spread_of(tab, n)
            d, tot = span_dim(S, n)
            st[(d, tot - d, selfadjoint_dim(tab, n), closed_pairs(S))] += 1
            pos = set(x for x in m if x > 0)
            s.add_clause([-v if v in pos else v for v in free])
    print("n=6 samples (span, dim Sigma^perp, self-adjoint dim, closed pairs): count")
    for k, v in sorted(st.items()):
        print(" ", k, v)
    n = 8
    D = list({frozenset({0, x, omega(x, n), x ^ omega(x, n)}) for x in range(1, 1 << n)})
    d, tot = span_dim(D, n)
    print(f"Desarguesian spread n=8: lines={len(D)} span={d} dim Sigma^perp={tot - d}")


if __name__ == "__main__":
    main()
