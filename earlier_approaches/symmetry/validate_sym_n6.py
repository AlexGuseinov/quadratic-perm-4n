"""Validation of the symmetric-case pipeline at n = 6 (where solutions exist).

For every one of the 51 symmetry pairs at n = 6 it compares
  (1) the plain verdict (logs/sym_n6.jsonl: kissat on the equivariant space, no lemmas),
  (2) the verdict with the lemma checks and normal-form rows of sym_search.normal_form
      (CaDiCaL through PySAT),
and, for the involution pairs, (3) the alternatives of sym_p2.py (generalised to n = 6:
x = e_{2k}), including the output normalisation (A9) in its n = 6 form (k = kB = 2).
A pair passes if all decided verdicts agree; SAT answers are checked to be permutations
with differential uniformity 4.
Usage: python3 validate_sym_n6.py   (writes logs/validate_sym_n6.log)
"""
import json
import os
from itertools import combinations

from pysat.solvers import Solver

import os, sys  # repository layout: shared modules are in ../../core
sys.path.insert(1, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "core"))
from du4perm_sat import extract, lut, check
from sym_search import build_sym, cases, fixed_dim, normal_form

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")   # repository root (logs/ is there)
n = 6
pairs = list(combinations(range(n), 2))
pidx = {pq: t for t, pq in enumerate(pairs)}


def zero(i, j):
    return [(1 << (pidx[(min(i, j), max(i, j))] * n + t), 0) for t in range(n)]


def setv(i, j, target):
    return [(1 << (pidx[(min(i, j), max(i, j))] * n + t), target >> t & 1) for t in range(n)]


def solve(A, B, rows, budget=5_000_000):
    try:
        enc, c, ell, dim = build_sym(n, A, B, True, False, rows)
    except ValueError:
        return "UNSAT"
    with Solver(name="cadical195", bootstrap_with=enc.clauses) as s:
        s.conf_budget(budget)
        res = s.solve_limited()
        if res is None:
            return "UNKNOWN"
        if res:
            C, L = extract(s.get_model(), c, ell, n)
            ok, du = check(lut(C, L, n), n)
            assert ok and du == 4, "SAT answer is not a DU-4 permutation"
            return "SAT"
        return "UNSAT"


plain = {}
for line in open(os.path.join(ROOT, "logs", "sym_n6.jsonl")):
    r = json.loads(line)
    if r.get("perm", True) and not r.get("nf") and not r.get("reduce"):
        plain[r["case"]] = r["status"]

out = []
fails = 0
for idx, (p, la, A, lb, B) in enumerate(cases(n)):
    verdicts = {"plain": plain.get(idx, "MISSING")}
    ok, why, rows = normal_form(p, A, B, n, la, lb)
    verdicts["normal_form"] = "UNSAT" if not ok else solve(A, B, rows)
    if p == 2 and ok:
        k, kB = n - fixed_dim(A, n), n - fixed_dim(B, n)
        if 2 * k < n:
            x = 2 * k
            alts = [zero(x, 0), zero(x, x + 1)]
            res = [solve(A, B, rows + a) for a in alts]
            verdicts["alternatives"] = "SAT" if "SAT" in res else ("UNSAT" if all(r == "UNSAT" for r in res) else "UNKNOWN")
            if k == kB == 2:      # (A9) at n = 6: w_i = beta(e4, e_{2i+1})
                outn = [zero(4, 0) + setv(4, 3, 1 << 1),
                        zero(4, 5) + setv(4, 1, 1 << 1) + setv(4, 3, 1 << 3)]
                res = [solve(A, B, rows + a) for a in outn]
                verdicts["output_norm"] = "SAT" if "SAT" in res else ("UNSAT" if all(r == "UNSAT" for r in res) else "UNKNOWN")
        else:
            alts = [zero(0, 2), zero(0, 1)]
            res = [solve(A, B, rows + a) for a in alts]
            verdicts["alternatives"] = "SAT" if "SAT" in res else ("UNSAT" if all(r == "UNSAT" for r in res) else "UNKNOWN")
    decided = {v for v in verdicts.values() if v in ("SAT", "UNSAT")}
    status = "PASS" if len(decided) == 1 and "UNKNOWN" not in verdicts.values() and "TIMEOUT" not in verdicts.values() else (
        "FAIL" if len(decided) > 1 else "UNDECIDED")
    fails += status == "FAIL"
    line = f"{idx:3d} p={p} A={la:<10s} B={lb:<10s} {status:9s} {verdicts} {why}"
    out.append(line)
    print(line, flush=True)
summary = f"pairs: {len(out)}, FAIL: {fails}, UNDECIDED: {sum('UNDECIDED' in l for l in out)}"
print(summary)
with open(os.path.join(ROOT, "logs", "validate_sym_n6.log"), "w") as f:
    f.write("\n".join(out + [summary]) + "\n")
