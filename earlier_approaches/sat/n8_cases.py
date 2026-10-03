"""Case split for the general (non-symmetric) n = 8 problem by the local structure of the
kernel spread.

For kernel lines L, M let k(L, M) be the number of kernel lines inside L + M (4-dim).
  * k = 5 ("closed pair") is impossible for n = 8 (lemma in NOTES_lemmas.md, "No closed
    pair"): proved by hand, no computation needed.
  * k = 4 is impossible in general (quadratic-form argument in NOTES_lemmas.md).
So every pair has k in {2, 3}. Cases:
  B1: some pair has k = 3.  WLOG (L1, L2) = (<e0,e1>, <e2,e3>) and the third line is
      N = {e0+e2, e1+e3, e0+e1+e2+e3} (choose e2 = phi(e0), e3 = phi(e1) for the
      isomorphism phi: L1 -> L2 whose graph is N); SB1-SB4 stay valid (for SB4 replace
      e1 by e0+e1 and e3 by e2+e3 together; SB2 is restored inside <e'0..e'5>).
  B2: every pair has k = 2; in particular the 6 pairs of basis lines.
In both cases the 6 pairs of basis lines are not closed; in B2 no cross point x = l + m
(l in Li, m in Lj) has its kernel line inside Li + Lj.

Usage: python3 n8_cases.py n CASE out.cnf [--noperm] [--local]   (CASE = B1 or B2)
       --local: s_a >= 2 witnesses only for a in U = <e0..e3> (a relaxation)
       --d4: add "k = 2 => dim span beta(L, M) = 4" for the basis pairs (valid for n = 8)
       python3 n8_cases.py n CASE --solve [--noperm]      (small n, CaDiCaL)
"""
import sys
from itertools import combinations

import os, sys  # repository layout: shared modules are in ../../core
sys.path.insert(1, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "core"))
from du4perm_sat import build, extract, lut, check, T, F_


def beta_vec(enc, cij, n, x, y):
    """Literals of beta(x, y) for vectors x, y (bitmasks over the basis)."""
    xs = [i for i in range(n) if x >> i & 1]
    ys = [j for j in range(n) if y >> j & 1]
    return [enc.xor([cij(i, j, t) for i in xs for j in ys]) for t in range(n)]


def nonzero(enc, lits):
    """A literal that is true iff some literal in lits is true."""
    lits = [l for l in lits if l is not F_]
    if any(l is T for l in lits):
        return T
    if not lits:
        return F_
    o = enc.new()
    enc.clauses.append([-o] + lits)
    for l in lits:
        enc.clauses.append([o, -l])
    return o


def and_all(enc, lits):
    acc = T
    for l in lits:
        acc = enc.and2(acc, l)
    return acc


def rank3_literal(enc, cij, n, x, U):
    """True iff beta(x, .) restricted to the 4-space U (list of 4 basis vectors, x in U)
    has kernel exactly <x>, i.e. the kernel line of x is not inside U."""
    span = {0}
    for v in U:
        span |= {s ^ v for s in span}
    comp = []
    cur = {0, x}
    for v in U:
        if v not in cur:
            comp.append(v)
            cur |= {s ^ v for s in cur}
        if len(comp) == 3:
            break
    vecs = [beta_vec(enc, cij, n, x, u) for u in comp]
    conds = []
    for r in range(1, 4):
        for T_ in combinations(range(3), r):
            comb = [enc.xor([vecs[i][t] for i in T_]) for t in range(n)]
            conds.append(nonzero(enc, comb))
    return and_all(enc, conds)


def basis_lines(n):
    return [(1 << (2 * i), 1 << (2 * i + 1)) for i in range(n // 2)]


def cross_points(L, M):
    return [l ^ m for l in (L[0], L[1], L[0] ^ L[1]) for m in (M[0], M[1], M[0] ^ M[1])]


def independent4(enc, cij, n, L, M):
    """Literal: the four vectors beta(l_i, m_j) (l_i in L, m_j in M basis) are independent,
    i.e. span beta(L, M) has dimension 4."""
    vecs = [beta_vec(enc, cij, n, l, m) for l in L for m in M]
    conds = []
    for r in range(1, 5):
        for T_ in combinations(range(4), r):
            comb = [enc.xor([vecs[i][t] for i in T_]) for t in range(n)]
            conds.append(nonzero(enc, comb))
    return and_all(enc, conds)


def make_extra(case, d4=False):
    """d4: add the consequence of the global count (NOTES, 'pair types'): for n = 8 every
    pair of kernel lines with k = 2 has dim span beta(L, M) = 4."""
    def extra(enc, cij, n):
        lines = basis_lines(n)
        for a, b in combinations(range(len(lines)), 2):
            L, M = lines[a], lines[b]
            U = [L[0], L[1], M[0], M[1]]
            if case == "B2":
                for x in cross_points(L, M):                 # every K_x leaves U
                    enc.clause([rank3_literal(enc, cij, n, x, U)])
                if d4:
                    enc.clause([independent4(enc, cij, n, L, M)])
            else:                                            # at least: not closed
                r3 = [rank3_literal(enc, cij, n, x, U) for x in cross_points(L, M)]
                enc.clause(r3)
                if d4 and (a, b) != (0, 1):                  # k = 2 (all r3 true) => d = 4
                    enc.clause([enc.neg(and_all(enc, r3)), independent4(enc, cij, n, L, M)])
        if case == "B1":
            e = [1 << i for i in range(n)]
            for lit in beta_vec(enc, cij, n, e[0] ^ e[2], e[1] ^ e[3]):   # N kernel line
                enc.clause([enc.neg(lit)])
            # (K4): with omega: e0 -> e1 -> e0+e1, e2 -> e3 -> e2+e3 (L, M, N are F_4-points),
            # q(x) = beta(x, omega x) takes one common value v on the six points of the other
            # two F_4-points P = F_4(e0+e3), Q = F_4(e0+e2+e3).
            def om(x):
                y = 0
                for i in range(4):
                    if x >> i & 1:
                        y ^= (1 << (i + 1)) if i % 2 == 0 else ((1 << (i - 1)) | (1 << i))
                return y
            PQ = []
            for base in (e[0] ^ e[3], e[0] ^ e[2] ^ e[3]):
                PQ += [base, om(base), om(om(base))]
            ref = beta_vec(enc, cij, n, PQ[0], om(PQ[0]))
            for p in PQ[1:]:
                other = beta_vec(enc, cij, n, p, om(p))
                for t in range(n):
                    enc.clause([enc.neg(enc.xor2(ref[t], other[t]))])
    return extra


def main():
    n, case = int(sys.argv[1]), sys.argv[2]
    perm = "--noperm" not in sys.argv
    a_list = list(range(1, 16)) if "--local" in sys.argv else None
    enc, c, ell = build(n, sb=True, sb3=True, perm=perm, extra=make_extra(case, "--d4" in sys.argv),
                        a_list=a_list)
    head = f"n={n} case={case} perm={perm} vars={enc.pool.top} clauses={len(enc.clauses)}"
    if "--solve" in sys.argv:
        from pysat.solvers import Solver
        print(head, flush=True)
        with Solver(name="cadical195", bootstrap_with=enc.clauses) as s:
            if s.solve():
                C, L = extract(s.get_model(), c, ell, n)
                print("SAT", check(lut(C, L, n), n) if perm else "(no permutation condition)")
            else:
                print("UNSAT")
        return
    out = sys.argv[3]
    with open(out, "w") as f:
        f.write(f"c {head}\n")
        f.write(f"p cnf {enc.pool.top} {len(enc.clauses)}\n")
        for cl in enc.clauses:
            f.write(" ".join(map(str, cl)) + " 0\n")
    print(head)


if __name__ == "__main__":
    main()
