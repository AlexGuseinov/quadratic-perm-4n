"""Case p = 5, A = B = f0+f0 (fixed-point-free) for n = 8, split into alternatives.

V = F_16^2 with A = multiplication by zeta (order 5); C(A) = GL_{F_16}(V) is transitive
on nonzero vectors, so WLOG consider the kernel line K_{e0} = <e0, v>:
  (a) v in F_16 e0 = <e0, e1, e2, e3>: 7 lines through e0 in this 4-space;
  (b) v outside F_16 e0: the stabiliser of e0 in GL_{F_16}(V) is transitive on the
      vectors outside F_16 e0, so WLOG v = e4.
Output: C(B) = GL_{F_16}(W) is transitive on nonzero vectors; beta(e0, e_j0) != 0 for
e_j0 outside K_{e0}, so WLOG beta(e0, e_j0) = e0 (j0 = 1, or 2 if e1 is in K_{e0}).
Each alternative is enumerated exactly (beta_enum + linear-part SAT, sym_twostage).
"""
import json, time, os
from itertools import combinations
import os, sys  # repository layout: shared modules are in ../../core
sys.path.insert(1, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "core"))
from sym_search import cases, normal_form, affine_space
from sym_twostage import stage1, stage2
from du4perm_sat import lut, check
n = 8
p, la, A, lb, B = cases(n)[44]
pairs = list(combinations(range(n), 2)); pidx = {pq: t for t, pq in enumerate(pairs)}
def u(i, j, t): return pidx[(min(i, j), max(i, j))] * n + t
def eq(i, vec, target):
    out = []
    for t in range(n):
        m = 0
        for j in range(n):
            if vec >> j & 1 and j != i:
                m ^= 1 << u(i, j, t)
        out.append((m, target >> t & 1))
    return out
alts, seen = [], set()
for v in range(2, 16):
    key = frozenset({v, v ^ 1})
    if key not in seen:
        seen.add(key); alts.append(v)
alts.append(1 << 4)
ok, why, base = normal_form(p, A, B, n, la, lb)
assert ok
log = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "logs", "sym_case44.jsonl")
for v in alts:
    j0 = 2 if v in (2, 3) else 1
    sol = affine_space(A, B, n, base + eq(0, v, 0) + eq(0, 1 << j0, 1))
    t0 = time.time()
    if sol is None:
        rec = {"v": v, "status": "UNSAT", "note": "inconsistent"}
    else:
        offset, basis = sol[0], sol[1]
        betas = stage1(n, basis, offset)
        found = None
        for C in betas:
            L = stage2(n, C)
            if L is not None:
                found = check(lut(C, L, n), n); break
        rec = {"v": v, "j0": j0, "bits": len(basis), "betas": len(betas),
               "status": "SAT" if found else "UNSAT", "check": found}
    rec.update(seconds=round(time.time() - t0, 1), when=time.strftime("%Y-%m-%d %H:%M:%S"))
    with open(log, "a") as f:
        f.write(json.dumps(rec) + "\n")
    print(rec, flush=True)
