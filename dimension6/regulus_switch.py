"""Regulus switching of the Desarguesian line spread and the space Sigma^perp.

For a line spread Sigma of PG(n-1, 2), Sigma^perp = {alternating B : B(L, L) = 0
for all L in Sigma}. A quadratic DU-4 permutation needs Sigma^perp to contain an
n-dimensional subspace of constant rank n-2 (L2, L3).

Experiment: start from the Desarguesian spread D (F_4-points), switch one or more
reguli (3 lines of D inside a closed 4-space -> the opposite regulus), and report
dim Sigma^perp, the rank distribution of Sigma^perp, and the number of closed pairs.
Usage: python3 regulus_switch.py n [switches]
"""

import random
import sys
from collections import Counter
from itertools import combinations

import os, sys  # repository layout: shared modules are in ../core
sys.path.insert(1, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "core"))
from quadtools import rank_rows


def omega(x, n):
    y = 0
    for t in range(0, n, 2):
        a, b = (x >> t) & 1, (x >> (t + 1)) & 1
        y |= (b << t) | ((a ^ b) << (t + 1))
    return y


def desarguesian(n):
    return list({frozenset({x, omega(x, n), x ^ omega(x, n)}) for x in range(1, 1 << n)})


def span(pts_sets):
    s = {0}
    for L in pts_sets:
        for p in L:
            s |= {p ^ x for x in s}
    return s


def opposite_regulus(R):
    """R: three pairwise skew lines in a 4-space U. Return the 3 transversal lines."""
    U = span(R)
    pts = [p for p in U if p]
    covered = set().union(*R)
    trans = set()
    for a, b in combinations(pts, 2):
        L = frozenset({a, b, a ^ b})
        if L <= covered and all(len(L & M) == 1 for M in R):
            trans.add(L)
    return list(trans)


def perp_space(S, n):
    """Basis of Sigma^perp as bitmasks over pairs (i<j)."""
    idx = list(combinations(range(n), 2))
    rows = []
    for L in S:
        u, v = sorted(L)[:2]
        m = 0
        for k, (i, j) in enumerate(idx):
            if ((u >> i) & 1) & ((v >> j) & 1) ^ ((u >> j) & 1) & ((v >> i) & 1):
                m |= 1 << k
        rows.append(m)
    P = len(idx)
    # null space of rows (as constraints <row, B> = 0)
    piv = {}
    for r in rows:
        for c, pr in piv.items():
            if r >> c & 1:
                r ^= pr
        if r:
            c = (r & -r).bit_length() - 1
            for c2 in list(piv):
                if piv[c2] >> c & 1:
                    piv[c2] ^= r
            piv[c] = r
    free = [c for c in range(P) if c not in piv]
    basis = []
    for f in free:
        v = 1 << f
        for c, pr in piv.items():
            if pr >> f & 1:
                v |= 1 << c
        basis.append(v)
    return basis, idx


def form_rank(mask, n, idx):
    rows = [0] * n
    for k, (i, j) in enumerate(idx):
        if mask >> k & 1:
            rows[i] |= 1 << j
            rows[j] |= 1 << i
    return rank_rows(rows)


def rank_distribution(basis, n, idx):
    d = len(basis)
    cnt = Counter()
    cur = 0
    for g in range(1, 1 << d):
        cur ^= basis[(g & -g).bit_length() - 1]
        cnt[form_rank(cur, n, idx)] += 1
    return dict(sorted(cnt.items()))


def closed_pairs(S):
    pt = {p: i for i, L in enumerate(S) for p in L}
    c = 0
    for A, B in combinations(S, 2):
        U = span([A, B])
        if all(S[pt[p]] <= U for p in U if p):
            c += 1
    return c


def closed_reguli(S):
    """All reguli R (3 lines of S) contained in a closed 4-space spanned by two lines."""
    pt = {p: i for i, L in enumerate(S) for p in L}
    out = set()
    for A, B in combinations(S, 2):
        U = span([A, B])
        inside = {pt[p] for p in U if p}
        if all(S[i] <= U for i in inside) and len(inside) == 5:
            for tri in combinations(sorted(inside), 3):
                out.add(tri)
    return [tuple(S[i] for i in tri) for tri in out]


def main():
    n = int(sys.argv[1])
    switches = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    trials = int(sys.argv[3]) if len(sys.argv) > 3 else 20
    D = desarguesian(n)
    b, idx = perp_space(D, n)
    print(f"n={n} Desarguesian: lines={len(D)} dim Sigma^perp={len(b)} ranks={rank_distribution(b, n, idx) if len(b) <= 14 else 'skipped'}")
    stats = Counter()
    random.seed(1)
    for t in range(trials):
        S = list(D)
        for _ in range(switches):
            regs = closed_reguli(S)
            if not regs:
                break
            R = random.choice(regs)
            S = [L for L in S if L not in R] + opposite_regulus(list(R))
        pts = [p for L in S for p in L]
        assert len(pts) == (1 << n) - 1 and len(set(pts)) == len(pts)
        b, idx = perp_space(S, n)
        rd = rank_distribution(b, n, idx) if len(b) <= 14 else "big"
        stats[(len(b), str(rd), closed_pairs(S))] += 1
    print(f"after {switches} random regulus switch(es), {trials} trials:")
    for k, v in sorted(stats.items()):
        print("  dim Sigma^perp, ranks, closed pairs:", k, "x", v)


if __name__ == "__main__":
    main()
