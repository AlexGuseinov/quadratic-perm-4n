"""Sample quadratic DU-4 permutations of F_2^6 from the exact model (du4perm_sat)
and test: are the component radicals exactly the lines of the isotropic spread?
Also records the Walsh/rank data and a simple EA-invariant (number of
distinct radicals) to see how many different structures occur."""

import random
import sys
from collections import Counter
from itertools import combinations

from pysat.solvers import Solver

import os, sys  # repository layout: shared modules are in ../core
sys.path.insert(1, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "core"))
from du4perm_sat import build, extract, lut, check
from quadtools import rank_rows


def analyse(table, n):
    N = 1 << n
    # polar rows of component b
    def rows(b):
        out = []
        for i in range(n):
            r = 0
            for j in range(n):
                if i != j:
                    v = table[(1 << i) ^ (1 << j)] ^ table[1 << i] ^ table[1 << j] ^ table[0]
                    if bin(b & v).count("1") & 1:
                        r |= 1 << j
            out.append(r)
        return out

    rads = {}
    for b in range(1, N):
        R = rows(b)
        assert rank_rows(R) == n - 2
        rads[b] = frozenset(a for a in range(N)
                            if all(bin(R[i] & a).count("1") % 2 == 0 for i in range(n)))
    # isotropic spread: ker beta(a, .)
    beta = lambda u, v: table[u ^ v] ^ table[u] ^ table[v] ^ table[0]
    spread = set()
    for a in range(1, N):
        spread.add(frozenset(v for v in range(N) if beta(a, v) == 0))
    return len(set(rads.values())), set(rads.values()) == spread


def main():
    n = 6
    samples = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    enc, c, ell = build(n)
    free = [v for p, vs in c.items() for v in vs if isinstance(v, int) and not isinstance(v, bool)]
    stats = Counter()
    with Solver(name="cadical195", bootstrap_with=enc.clauses) as s:
        for t in range(samples):
            random.seed(1000 + t)
            assum = [v if random.random() < 0.5 else -v for v in random.sample(free, 6)]
            if not s.solve(assumptions=assum) and not s.solve():
                break
            model = s.get_model()
            C, L = extract(model, c, ell, n)
            table = lut(C, L, n)
            perm, du = check(table, n)
            assert perm and du == 4
            stats[analyse(table, n)] += 1
            pos = set(x for x in model if x > 0)
            s.add_clause([-v if v in pos else v for v in free])
    print("(number of distinct radicals, radicals == isotropic spread): count")
    for k, v in sorted(stats.items()):
        print(" ", k, v)


if __name__ == "__main__":
    main()
