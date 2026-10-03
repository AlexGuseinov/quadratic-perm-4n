"""SAT search for constant-rank subspaces of alternating matrices over F_2.

Question: is there a d-dimensional subspace W of Alt(n, F_2) in which every
nonzero element has rank exactly r?  For (n, r, d) = (8, 6, 8) this is the
necessary condition for a quadratic 8-bit permutation with differential
uniformity 4 (hypothesis H says: no).

Encoding
  * basis matrices B_1..B_d, one Boolean variable per upper-triangle entry;
    B_1 is fixed to the canonical rank-r form e0^e1 + e2^e3 + ... (WLOG under
    congruence, since all rank-r alternating forms are equivalent);
  * for every nonzero c in F_2^d the entries of M_c = sum c_k B_k are XORs of
    basis entries (built incrementally along c);
  * rank(M) <= r  <=>  every (r+2)x(r+2) principal Pfaffian vanishes;
    rank(M) >= r  <=>  some r x r principal Pfaffian is nonzero.
    Pfaffians are expanded along the first index and memoised:
    Pf(S) = XOR_j m[i0][j] AND Pf(S \\ {i0, j}).

Usage: python3 constrank_sat.py n r d [time_limit_s]
"""

import sys
import time
from itertools import combinations

from pysat.formula import IDPool
from pysat.solvers import Solver


class Encoder:
    def __init__(self):
        self.pool = IDPool()
        self.clauses = []
        self.TRUE = self.pool.id("TRUE")
        self.clauses.append([self.TRUE])

    def new(self):
        return self.pool.id(("aux", len(self.pool.obj2id)))

    def and2(self, a, b):
        if a == self.TRUE:
            return b
        if b == self.TRUE:
            return a
        o = self.new()
        self.clauses += [[-o, a], [-o, b], [o, -a, -b]]
        return o

    def xor2(self, a, b):
        o = self.new()
        self.clauses += [[-o, a, b], [-o, -a, -b], [o, -a, b], [o, a, -b]]
        return o

    def xor(self, lits):
        acc = lits[0]
        for x in lits[1:]:
            acc = self.xor2(acc, x)
        return acc


def build(n, r, d, du4=False):
    """du4=True adds, for every a != 0, two distinct nonzero witnesses
    c with M_c a = 0 (s_a >= 2).  With constant rank n-2 and d = n the
    count identity sum_a (2^{s_a}-1) = 3(2^n-1) then forces s_a = 2 for all a,
    i.e. differential uniformity exactly 4 (necessary for a quadratic DU-4
    permutation)."""
    enc = Encoder()
    pairs = list(combinations(range(n), 2))
    pidx = {p: k for k, p in enumerate(pairs)}

    # basis variables; B_1 fixed canonical
    basis = []
    canon = {(2 * i, 2 * i + 1) for i in range(r // 2)}
    for k in range(d):
        row = []
        for p in pairs:
            v = enc.pool.id(("b", k, p))
            if k == 0:
                enc.clauses.append([v] if p in canon else [-v])
            row.append(v)
        basis.append(row)

    # entries of every combination M_c
    entries = {}
    for c in range(1, 1 << d):
        low = (c & -c).bit_length() - 1
        rest = c ^ (1 << low)
        if rest == 0:
            entries[c] = basis[low]
        else:
            entries[c] = [enc.xor2(entries[rest][q], basis[low][q])
                          for q in range(len(pairs))]

    for c in range(1, 1 << d):
        m = entries[c]
        memo = {}

        def pf(S):
            if not S:
                return enc.TRUE
            if len(S) == 2:
                return m[pidx[S]]
            if S in memo:
                return memo[S]
            i0 = S[0]
            terms = []
            for j in S[1:]:
                sub = tuple(x for x in S if x not in (i0, j))
                terms.append(enc.and2(m[pidx[(i0, j)]], pf(sub)))
            out = enc.xor(terms)
            memo[S] = out
            return out

        # rank <= r
        if r + 2 <= n:
            for S in combinations(range(n), r + 2):
                enc.clauses.append([-pf(S)])
        # rank >= r
        enc.clauses.append([pf(S) for S in combinations(range(n), r)])

    if du4:
        # (B_k a)_i = XOR_{j in a} b_k[i][j]
        def entry(k, i, j):
            if i == j:
                return None
            p = (i, j) if i < j else (j, i)
            return basis[k][pidx[p]]

        for a in range(1, 1 << n):
            Ba = []  # Ba[k][i]
            for k in range(d):
                col = []
                for i in range(n):
                    lits = [entry(k, i, j) for j in range(n)
                            if (a >> j) & 1 and j != i]
                    col.append(enc.xor(lits) if lits else None)
                Ba.append(col)
            wits = []
            for w in range(2):
                u = [enc.pool.id(("w", a, w, k)) for k in range(d)]
                wits.append(u)
                enc.clauses.append(u[:])  # u != 0
                for i in range(n):
                    terms = [enc.and2(u[k], Ba[k][i]) for k in range(d)
                             if Ba[k][i] is not None]
                    if terms:
                        enc.clauses.append([-enc.xor(terms)])
            # u != v
            diff = [enc.xor2(wits[0][k], wits[1][k]) for k in range(d)]
            enc.clauses.append(diff)
    return enc, basis, pairs


def decode(model, basis, pairs, n):
    pos = set(x for x in model if x > 0)
    mats = []
    for row in basis:
        rows = [0] * n
        for v, (i, j) in zip(row, pairs):
            if v in pos:
                rows[i] |= 1 << j
                rows[j] |= 1 << i
        mats.append(rows)
    return mats


def rank_rows(rows):
    rows = [x for x in rows if x]
    rk = 0
    while rows:
        p = rows.pop()
        rk += 1
        low = p & -p
        rows = [x ^ p if x & low else x for x in rows]
        rows = [x for x in rows if x]
    return rk


def combo(mats, c, n):
    acc = [0] * n
    for k in range(len(mats)):
        if c >> k & 1:
            acc = [x ^ y for x, y in zip(acc, mats[k])]
    return acc


def verify(mats, n, r):
    d = len(mats)
    return all(rank_rows(combo(mats, c, n)) == r for c in range(1, 1 << d))


def s_profile(mats, n):
    """s_a = dim{c : M_c a = 0} for every a != 0, returned as a histogram."""
    d = len(mats)
    hist = {}
    for a in range(1, 1 << n):
        cnt = 0
        for c in range(1 << d):
            M = combo(mats, c, n)
            if all(bin(M[i] & a).count("1") % 2 == 0 for i in range(n)):
                cnt += 1
        s = cnt.bit_length() - 1
        hist[s] = hist.get(s, 0) + 1
    return hist


def main():
    args = [x for x in sys.argv[1:] if not x.startswith("--")]
    du4 = "--du4" in sys.argv
    n, r, d = map(int, args[:3])
    limit = float(args[3]) if len(args) > 3 else None
    t0 = time.time()
    enc, basis, pairs = build(n, r, d, du4=du4)
    print(f"n={n} r={r} d={d} du4={du4}: vars={enc.pool.top} "
          f"clauses={len(enc.clauses)} (build {time.time() - t0:.1f}s)",
          flush=True)
    with Solver(name="cadical195", bootstrap_with=enc.clauses) as s:
        t1 = time.time()
        if limit:
            import threading
            timer = threading.Timer(limit, s.interrupt)
            timer.start()
            res = s.solve_limited(expect_interrupt=True)
            timer.cancel()
        else:
            res = s.solve()
        dt = time.time() - t1
        if res is True:
            mats = decode(s.get_model(), basis, pairs, n)
            ok = verify(mats, n, r)
            print(f"SAT in {dt:.1f}s; independent rank check: {ok}; "
                  f"s_a histogram: {s_profile(mats, n)}")
            for k, rows in enumerate(mats):
                print(f"  B{k + 1}: " + " ".join(f"{x:0{n}b}"[::-1] for x in rows))
        elif res is False:
            print(f"UNSAT in {dt:.1f}s")
        else:
            print(f"UNKNOWN (time limit {limit}s reached after {dt:.1f}s)")


if __name__ == "__main__":
    main()
