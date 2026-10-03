"""Two-stage exact search for one symmetry case of sym_search.py.

Stage 1 (beta_enum): all beta in the equivariant space with s_a = 2 for all a != 0
and every component of rank n-2 (both conditions depend on beta only).
Stage 2 (SAT, 64 variables for n = 8): for each surviving beta, is there a linear
part L with Q(a) + L(a) not in Im beta(a, .) for all a != 0 (Q the canonical
quadratic map with polar beta)?  That is exactly: F = Q + L is a permutation.
Stage 1 also gives the number of symmetric beta without the permutation condition.
Usage: python3 sym_twostage.py n [--nf] CASE [CASE ...]
  --nf: add the normal-form conditions of sym_search.normal_form (odd p)
"""
import json
import os
import subprocess
import sys
import time
from itertools import combinations

from pysat.solvers import Solver

import os, sys  # repository layout: shared modules are in ../../core
sys.path.insert(1, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "core"))
from du4perm_sat import Enc, T, F_, lut, check
from sym_search import cases, equivariant_basis, normal_form, affine_space

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")   # repository root (logs/ is there)


def stage1(n, basis, offset=0):
    pairs = list(combinations(range(n), 2))
    lines = [f"{n} {len(basis)} 0"]
    for v in [offset] + basis:
        lines.append("".join(f"{(v >> (n * p)) & ((1 << n) - 1):02x}"
                             for p in range(len(pairs))))
    r = subprocess.run([os.path.join(HERE, "beta_enum")], input="\n".join(lines) + "\n",
                       capture_output=True, text=True, check=True)
    out = r.stdout.split("\n")
    count = int(out[-2].split()[1])
    betas = []
    for line in out[:-2]:
        vals = [int(line[2 * p:2 * p + 2], 16) for p in range(len(pairs))]
        betas.append(dict(zip(pairs, vals)))
    assert len(betas) == count
    return betas


def stage2(n, C):
    """Find L (list of F(e_i)) making F = Q + L a permutation, or None."""
    N = 1 << n

    def cij(i, j):
        return 0 if i == j else C[(min(i, j), max(i, j))]

    enc = Enc()
    ell = [[enc.var(("l", i, k)) for k in range(n)] for i in range(n)]
    for a in range(1, N):
        bits = [i for i in range(n) if a >> i & 1]
        rows = [0] * n                       # beta(a, e_j)
        for j in range(n):
            for i in bits:
                rows[j] ^= cij(i, j)
        W = [w for w in range(1, N) if all(bin(w & r).count("1") % 2 == 0 for r in rows)]
        assert len(W) == 3                   # (Im beta(a, .))^perp has dimension 2
        q = 0
        for i, j in combinations(bits, 2):
            q ^= cij(i, j)
        lits = []
        for w in W[:2]:
            const = bin(w & q).count("1") % 2
            x = enc.xor([ell[i][k] for i in bits for k in range(n) if w >> k & 1])
            lits.append(enc.neg(x) if const else x)
        enc.clause(lits)                     # <w, Q(a) + L(a)> = 1 for w1 or w2
    with Solver(name="cadical195", bootstrap_with=enc.clauses) as s:
        if not s.solve():
            return None
        pos = set(x for x in s.get_model() if x > 0)
        return [sum((1 << k) for k in range(n) if ell[i][k] in pos) for i in range(n)]


def main():
    n = int(sys.argv[1])
    cs = cases(n)
    nf = "--nf" in sys.argv
    out = os.path.join(ROOT, "logs", f"sym2_n{n}{'_nf' if nf else ''}.jsonl")
    for idx in [int(x) for x in sys.argv[2:] if not x.startswith("--")]:
        p, la, A, lb, B = cs[idx]
        offset = 0
        if nf:
            ok, why, rows = normal_form(p, A, B, n, la, lb)
            sol = affine_space(A, B, n, rows) if ok else None
            if sol is None:
                rec = {"n": n, "case": idx, "p": p, "A": la, "B": lb, "status": "UNSAT",
                       "reason": why or "normal-form equations inconsistent",
                       "when": time.strftime("%Y-%m-%d %H:%M:%S")}
                with open(out, "a") as f:
                    f.write(json.dumps(rec) + "\n")
                print(f"case {idx:3d} p={p:<3d} A={la:<12s} B={lb:<12s} UNSAT ({rec['reason']})", flush=True)
                continue
            offset, basis = sol[0], sol[1]
        else:
            basis = equivariant_basis(A, B, n)[0]
        t0 = time.time()
        betas = stage1(n, basis, offset)
        t1 = time.time()
        found = None
        for C in betas:
            L = stage2(n, C)
            if L is not None:
                table = lut(C, L, n)
                ok, du = check(table, n)
                found = {"lut": table, "check_perm": ok, "check_du": du}
                break
        rec = {"n": n, "case": idx, "p": p, "A": la, "B": lb, "free_beta_bits": len(basis),
               "stage1_betas": len(betas), "stage1_seconds": round(t1 - t0, 1),
               "stage2_seconds": round(time.time() - t1, 1),
               "status": "SAT" if found else "UNSAT", "when": time.strftime("%Y-%m-%d %H:%M:%S")}
        if found:
            rec.update(found)
        with open(out, "a") as f:
            f.write(json.dumps(rec) + "\n")
        print(f"case {idx:3d} p={p:<3d} A={la:<12s} B={lb:<12s} bits={len(basis):3d} "
              f"betas={len(betas):8d} {rec['status']} "
              f"({rec['stage1_seconds']}s + {rec['stage2_seconds']}s)"
              + (f" perm={found['check_perm']} DU={found['check_du']}" if found else ""),
              flush=True)


if __name__ == "__main__":
    main()
