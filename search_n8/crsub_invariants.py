"""GL-invariants of the subspaces W found by crsub_enum (small n), for validation.

For W (component space, n basis forms as pair masks) compute:
  * the kernel spread {K_a} (checks that every K_a is a line and that they form a spread);
  * the number of pairs of spread lines L, M with k(L, M) = number of spread lines in L + M,
    as a histogram {k: count};
  * the number of distinct radicals and how many radicals are spread lines;
  * whether a linear part L exists with Q + L a permutation (sym_twostage.stage2).
Usage: python3 crsub_invariants.py n solutions.txt [--max N] [--perm]
"""
import random
import sys
from collections import Counter
from itertools import combinations

import os, sys  # repository layout: sym_twostage is in ../earlier_approaches/symmetry
sys.path.insert(1, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "earlier_approaches", "symmetry"))
from sym_twostage import stage2


def pairs_of(n):
    return list(combinations(range(n), 2))


def form_row(mask, a, n, pairs):
    """w(a, .) as an n-bit row."""
    row = 0
    for p, (i, j) in enumerate(pairs):
        if mask >> p & 1:
            if a >> i & 1:
                row ^= 1 << j
            if a >> j & 1:
                row ^= 1 << i
    return row


def invariants(basis, n, perm=False):
    pairs = pairs_of(n)
    N = 1 << n
    # beta(a, x)_k = w_k(a, x)
    rows = {a: [form_row(w, a, n, pairs) for w in basis] for a in range(1, N)}
    spread = set()
    for a in range(1, N):
        K = frozenset(x for x in range(1, N) if all(bin(r & x).count("1") % 2 == 0 for r in rows[a]))
        assert len(K) == 3 and a in K, "kernel is not a line"
        spread.add(K)
    assert len(spread) == (N - 1) // 3, "kernel lines do not form a spread"
    lines = sorted(tuple(sorted(L)) for L in spread)
    hist = Counter()
    for L, M in combinations(lines, 2):
        U = {0}
        for v in L + M:
            U |= {u ^ v for u in U}
        k = sum(1 for K in lines if set(K) <= U)
        hist[k] += 1
    # radicals of all nonzero elements of W
    elems = [0]
    for w in basis:
        elems += [e ^ w for e in elems]
    rads = Counter()
    for e in elems[1:]:
        R = frozenset(x for x in range(1, N) if form_row(e, x, n, pairs) == 0)
        assert len(R) == 3
        rads[R] += 1
    rad_spread = sum(c for R, c in rads.items() if R in spread)
    out = {"k_hist": dict(sorted(hist.items())), "distinct_radicals": len(rads),
           "radicals_on_spread": rad_spread}
    if perm:
        C = {}
        for p, (i, j) in enumerate(pairs):
            C[(i, j)] = sum(((w >> p) & 1) << k for k, w in enumerate(basis))
        out["permutation"] = stage2(n, C) is not None
    return out


def main():
    n = int(sys.argv[1])
    lines = [l.split()[1:] for l in open(sys.argv[2]) if l.startswith("SOLUTION")]
    seen, sols = set(), []
    for l in lines:
        b = tuple(sorted((int(x) for x in l), reverse=True))
        if b not in seen:
            seen.add(b)
            sols.append(list(b))
    mx = int(sys.argv[sys.argv.index("--max") + 1]) if "--max" in sys.argv else len(sols)
    random.seed(1)
    if len(sols) > mx:
        sols = random.sample(sols, mx)
    perm = "--perm" in sys.argv
    tally = Counter()
    for b in sols:
        inv = invariants(b, n, perm)
        tally[repr(inv)] += 1
    print(f"distinct solutions: {len(seen)}, analysed: {len(sols)}")
    for k, v in sorted(tally.items(), key=lambda t: -t[1]):
        print(f"{v:6d}  {k}")


if __name__ == "__main__":
    main()
