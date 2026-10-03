"""Numerical check of the 'no closed pair' lemma ingredients on n = 6 solutions.

Lemma (n = 8). If two kernel lines L, M span a 4-space U that is a union of five kernel
lines, then: (i) every B_b restricted to U has rank 0 or 4 (U carries the F_4-structure of
its regular spread and beta(x, omega x) = 0 on U); (ii) the counting identity on U gives
dim S = 2, S = span beta(U, U); (iii) for b != 0 orthogonal to S, U is totally isotropic for
B_b, a form of rank n - 2, so dim(U & R_b) >= 4 - (n-2)/2; (iv) every a != 0 lies in exactly
3 radicals. For n = 8, (iii) gives |R_b & U| >= 1 for 63 values of b, while (iv) bounds the
number of pairs (a in U, b) with a in R_b by 45: contradiction. For n = 6, (iii) says R_b lies
in U for the 15 such b, i.e. exactly 45 pairs -- consistent, and closed pairs do occur.
This script checks (i)-(iv) and the n = 6 equality on sampled n = 6 solutions.
"""
import random
from itertools import combinations
from pysat.solvers import Solver
import os, sys  # repository layout: shared modules are in ../core
sys.path.insert(1, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "core"))
from du4perm_sat import build, extract, lut

n = 6
N = 1 << n


def span(vs):
    s = {0}
    for v in vs:
        s |= {x ^ v for x in s}
    return s


def rank_rows(rows):
    rows = list(rows); r = 0
    for bit in range(n):
        piv = next((x for x in rows if x >> bit & 1), None)
        if piv is None:
            continue
        rows = [x ^ piv if x >> bit & 1 else x for x in rows if x is not piv]
        r += 1
    return r


enc, c, ell = build(n, sb=True, sb3=True)
free = [v for p, vs in c.items() for v in vs if isinstance(v, int) and not isinstance(v, bool)]
checked = 0
with Solver(name="cadical195", bootstrap_with=enc.clauses) as s:
    for t in range(10):
        random.seed(500 + t)
        assum = [v if random.random() < 0.5 else -v for v in random.sample(free, 8)]
        if not s.solve(assumptions=assum):
            continue
        C, L = extract(s.get_model(), c, ell, n)
        tab = lut(C, L, n)
        beta = lambda u, v: tab[u ^ v] ^ tab[u] ^ tab[v]
        K = {a: frozenset(v for v in range(N) if beta(a, v) == 0) for a in range(1, N)}
        lines = list(set(K.values()))
        R = {}
        for b in range(1, N):
            R[b] = frozenset(a for a in range(N)
                             if all(bin(b & beta(a, v)).count("1") % 2 == 0 for v in range(N)))
        for Lx, Mx in combinations(lines, 2):
            U = span(list(Lx) + list(Mx))
            if sum(1 for X in lines if X <= U) != 5:
                continue
            S = span([beta(u, v) for u in U for v in U])
            assert len(S) == 4                                     # (ii) dim S = 2
            basis = [min(Lx - {0}), max(Lx - {0}), min(Mx - {0}), max(Mx - {0})]
            for b in range(1, N):
                Mb = [sum(((bin(b & beta(u, w)).count("1") & 1) << j) for j, w in enumerate(basis)) for u in basis]
                assert rank_rows(Mb) in (0, 4)                     # (i)
            perp = [b for b in range(1, N) if all(bin(b & x).count("1") % 2 == 0 for x in S)]
            assert len(perp) == 2 ** (n - 2) - 1
            pairs_count = sum(len((R[b] & U) - {0}) for b in perp)
            assert all(R[b] <= U for b in perp)                    # (iii) at n = 6
            assert pairs_count == 45                               # (iv) equality at n = 6
            checked += 1
print("closed pairs checked:", checked, "- all assertions hold")
