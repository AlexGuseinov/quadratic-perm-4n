"""Exact SAT model for: does F_2^n admit a quadratic permutation with
differential uniformity 4?  (Question D(n); main target n = 8.)

F(x) = sum_i x_i l_i + sum_{i<j} x_i x_j c_ij,   l_i, c_ij in F_2^n,  F(0) = 0.
beta(e_i, e_j) = c_ij; component b has alternating matrix B_b[i][j] = <b, c_ij>.

Constraints (each justified in NOTES_lemmas.md):
  * every component has rank exactly n-2                      (L2)
  * for every a != 0 there are two distinct nonzero b in S_a  (s_a >= 2, L2)
    (w1_a, w2_a witnesses: sum_k w_k B_k a = 0)
  * F is a permutation: for every a != 0, <w1_a, F(a)> = 1 or <w2_a, F(a)> = 1
    (D_a F(x) = beta(a, x) + F(a) is never 0 iff F(a) is not in Im beta(a, .)
     = S_a^perp; with s_a = 2, S_a = span(w1_a, w2_a))
Symmetry breaking (WLOG, proofs in NOTES_lemmas.md, "SB1/SB2"):
  SB1  the lines <e_{2i}, e_{2i+1}> belong to the isotropic spread:
       c_{2i,2i+1} = 0 for all i;
  SB2  beta(e_0, e_j) = u_{j-2} (unit vectors) for j = 2..n-1.

Usage: python3 du4perm_sat.py n [--sb3] [--noperm] [--desarg] [--closed] [--nosb] [--conflicts C]
"""

import sys
import time
from itertools import combinations

from pysat.formula import IDPool
from pysat.solvers import Solver

T, F_ = True, False  # constant literals


class Enc:
    def __init__(self):
        self.pool = IDPool()
        self.clauses = []

    def var(self, name):
        return self.pool.id(name)

    def new(self):
        return self.pool.id(("aux", self.pool.top + 1))

    def neg(self, a):
        if a is T:
            return F_
        if a is F_:
            return T
        return -a

    def and2(self, a, b):
        if a is F_ or b is F_:
            return F_
        if a is T:
            return b
        if b is T:
            return a
        o = self.new()
        self.clauses += [[-o, a], [-o, b], [o, -a, -b]]
        return o

    def xor2(self, a, b):
        if a is F_:
            return b
        if b is F_:
            return a
        if a is T:
            return self.neg(b)
        if b is T:
            return self.neg(a)
        o = self.new()
        self.clauses += [[-o, a, b], [-o, -a, -b], [o, -a, b], [o, a, -b]]
        return o

    def xor(self, lits):
        acc = F_
        for x in lits:
            acc = self.xor2(acc, x)
        return acc

    def clause(self, lits):
        """Add a clause; constants are folded (True satisfies, False dropped)."""
        out = []
        for x in lits:
            if x is T:
                return
            if x is F_:
                continue
            out.append(x)
        if not out:
            raise ValueError("empty clause: instance trivially UNSAT")
        self.clauses.append(out)


def build(n, sb=True, sb3=False, perm=True, desarg=False, closed=False, extra=None,
          a_list=None):
    """sb3: additionally F(e_0) = u_{n-2} (SB3) and bit n-2 of F(e_1) = 0 (SB4).
    closed: case (A) of the n = 8 split -- U = <e0, e1, e2, e3> is a union of five
    kernel lines, normalised to the F_4-spread of omega (e0 -> e1, e2 -> e3), i.e.
    beta(x, omega x) = 0 for x in U (the Desarguesian constraints restricted to U)."""
    enc = Enc()
    pairs = list(combinations(range(n), 2))
    c = {}
    for (i, j) in pairs:
        c[(i, j)] = [enc.var(("c", i, j, k)) for k in range(n)]
    if sb:
        for i in range(0, n, 2):
            c[(i, i + 1)] = [F_] * n                                  # SB1
        for j in range(2, n):
            c[(0, j)] = [T if k == j - 2 else F_ for k in range(n)]   # SB2
    ell = [[enc.var(("l", i, k)) for k in range(n)] for i in range(n)]
    if sb and sb3:
        ell[0] = [T if k == n - 2 else F_ for k in range(n)]          # SB3
        ell[1][n - 2] = F_                                             # SB4

    def cij(i, j, k):
        if i == j:
            return F_
        return c[(i, j) if i < j else (j, i)][k]

    ddim = n if desarg else (4 if closed else 0)
    if ddim:
        # Desarguesian isotropic spread: the lines <e_{2t}, e_{2t+1}> are F_4 e_{2t}
        # with omega e_{2t} = e_{2t+1}, omega e_{2t+1} = e_{2t} + e_{2t+1};
        # require beta(x, omega x) = 0 for all x (linear conditions on beta).
        def Me(k):  # omega e_k as a set of basis indices
            return [k + 1] if k % 2 == 0 else [k - 1, k]
        def beta_e_M(i, k, t):  # t-th coordinate of beta(e_i, omega e_k)
            return enc.xor([cij(i, j, t) for j in Me(k)])
        for t in range(n):
            for i in range(ddim):
                enc.clause([enc.neg(beta_e_M(i, i, t))])
                for k in range(i + 1, ddim):
                    enc.clause([enc.neg(enc.xor2(beta_e_M(i, k, t), beta_e_M(k, i, t)))])

    if extra is not None:
        extra(enc, cij, n)          # additional case constraints (n8_cases.py)
    add_core(enc, n, c, ell, perm, a_list=a_list)
    return enc, c, ell


def add_core(enc, n, c, ell, perm=True, a_list=None, b_list=None):
    """Rank-(n-2), s_a >= 2 witness and permutation constraints for given
    literals c[(i, j)][k] (beta(e_i, e_j)) and ell[i][k] (F(e_i)).
    a_list / b_list: restrict the witness / rank constraints to these a / b
    (used when a symmetry makes the other constraints redundant); default all."""
    pairs = list(combinations(range(n), 2))

    def cij(i, j, k):
        if i == j:
            return F_
        return c[(i, j) if i < j else (j, i)][k]

    # --- rank exactly n-2 for every component b
    pidx = {p: t for t, p in enumerate(pairs)}
    for b in (b_list if b_list is not None else range(1, 1 << n)):
        m = [enc.xor([c[p][k] for k in range(n) if b >> k & 1]) for p in pairs]
        memo = {}

        def pf(S):
            if not S:
                return T
            if len(S) == 2:
                return m[pidx[S]]
            if S in memo:
                return memo[S]
            i0 = S[0]
            terms = []
            for j in S[1:]:
                sub = tuple(x for x in S if x not in (i0, j))
                terms.append(enc.and2(m[pidx[(i0, j)]], pf(sub)))
            memo[S] = enc.xor(terms)
            return memo[S]

        enc.clause([enc.neg(pf(tuple(range(n))))])                      # rank <= n-2
        enc.clause([pf(S) for S in combinations(range(n), n - 2)])     # rank >= n-2

    # --- witnesses and permutation condition
    wit = {}
    for a in (a_list if a_list is not None else range(1, 1 << n)):
        bits = [i for i in range(n) if a >> i & 1]
        Ba = [[enc.xor([cij(i, j, k) for j in bits]) for i in range(n)]
              for k in range(n)]                 # Ba[k][i] = (B_k a)_i
        Fa = [enc.xor([ell[i][k] for i in bits] +
                      [cij(i, j, k) for i, j in combinations(bits, 2)])
              for k in range(n)]                 # F(a)_k
        ws, ps = [], []
        for w in range(2):
            u = [enc.var(("w", a, w, k)) for k in range(n)]
            ws.append(u)
            enc.clause(u[:])                                            # u != 0
            for i in range(n):
                enc.clause([enc.neg(enc.xor([enc.and2(u[k], Ba[k][i])
                                             for k in range(n)]))])     # M_u a = 0
            ps.append(enc.xor([enc.and2(u[k], Fa[k]) for k in range(n)]))
        enc.clause([enc.xor2(ws[0][k], ws[1][k]) for k in range(n)])  # u != v
        if perm:
            enc.clause(ps)                                              # permutation
        wit[a] = ws
    return wit


def extract(model, c, ell, n):
    pos = set(x for x in model if x > 0)

    def val(lit):
        if lit is T:
            return 1
        if lit is F_:
            return 0
        return 1 if lit in pos else 0

    C = {p: sum(val(v) << k for k, v in enumerate(vs)) for p, vs in c.items()}
    L = [sum(val(v) << k for k, v in enumerate(row)) for row in ell]
    return C, L


def lut(C, L, n):
    table = []
    for x in range(1 << n):
        y = 0
        for i in range(n):
            if x >> i & 1:
                y ^= L[i]
                for j in range(i + 1, n):
                    if x >> j & 1:
                        y ^= C[(i, j)]
        table.append(y)
    return table


def check(table, n):
    N = 1 << n
    perm = len(set(table)) == N
    du = 0
    for a in range(1, N):
        cnt = [0] * N
        for x in range(N):
            cnt[table[x ^ a] ^ table[x]] += 1
        du = max(du, max(cnt))
    return perm, du


def main():
    n = int(sys.argv[1])
    seed = int(sys.argv[sys.argv.index("--seed") + 1]) if "--seed" in sys.argv else 0
    budget = int(sys.argv[sys.argv.index("--conflicts") + 1]) if "--conflicts" in sys.argv else 0
    nosb = "--nosb" in sys.argv
    sb3 = "--sb3" in sys.argv
    perm = "--noperm" not in sys.argv
    desarg = "--desarg" in sys.argv
    closed = "--closed" in sys.argv
    t0 = time.time()
    enc, c, ell = build(n, sb=not nosb, sb3=sb3, perm=perm, desarg=desarg, closed=closed)
    print(f"n={n} sb={not nosb} sb3={sb3} perm={perm} desarg={desarg}: "
          f"vars={enc.pool.top} clauses={len(enc.clauses)} "
          f"(build {time.time() - t0:.1f}s)", flush=True)
    with Solver(name="cadical195", bootstrap_with=enc.clauses) as s:
        t1 = time.time()
        if budget:
            s.conf_budget(budget)
            res = s.solve_limited()
        else:
            res = s.solve()
        dt = time.time() - t1
        if res is True:
            C, L = extract(s.get_model(), c, ell, n)
            table = lut(C, L, n)
            perm, du = check(table, n)
            print(f"SAT in {dt:.1f}s: independent check permutation={perm} DU={du}")
            print("LUT:", table)
        elif res is False:
            print(f"UNSAT in {dt:.1f}s")
        else:
            print(f"UNKNOWN after {dt:.1f}s (conflict budget {budget})")


if __name__ == "__main__":
    main()
