"""Remaining involution cases for n = 8 (p = 2), split into alternatives, each solved
with kissat (DRAT-checked when short). A = J2^k (blocks on (e0,e1),(e2,e3),...), so
V1 = ker(A - I) and U_A = Im(A - I) = <e0, e2, ..., e_{2k-2}>.
  k = 3: V1 = U_A + <e6, e7>. C(A) is transitive on V1 \\ U_A (it contains
         e6 -> e6 + u, u in U_A, and GL(<e6, e7>)), so WLOG x = e6; K_{e6} lies in V1
         (N10) and is <e6, u> (u in U_A \\ 0; WLOG u = e0 via GL(3,2) lifted to C(A)) or
         <e6, e7 + u> (WLOG u = 0 via e7 -> e7 + u). Alternatives beta(e6, e0) = 0 and
         beta(e6, e7) = 0.
  k = 4: V1 = U_A = <e0, e2, e4, e6>; C(A) = GL(4, F_2[t]/t^2) acts on V1 as GL(4, 2), so
         WLOG x = e0 and K_{e0} = <e0, e2> (inside V1) or <e0, e1 + v> (v in V1; WLOG v = 0
         via I + (A - I)X in C(A)). Alternatives beta(e0, e2) = 0 and beta(e0, e1) = 0.
         If kB = 2, V1 is a union of kernel lines (N10), a regular spread of PG(3, 2);
         WLOG its F_4-structure is e0 -> e2 -> e0+e2, e4 -> e6 -> e4+e6.
  --outnorm (k = kB = 3, n = 8): also normalise the output. With eps = B - I, W is the
         F_2[eps]-module F_2[eps]^3 + F_2^2 and C(B) = Aut(W). The vectors
         w_i = beta(e6, e_{2i+1}) satisfy eps w_i = beta(e6, e_{2i}) (N10 relation). In the
         alternative K_{e6} = <e6, e7>, eps w_0, eps w_1, eps w_2 are independent, so the
         w_i span a pure free submodule (a direct summand): WLOG w_i = e'_{2i+1}. In the
         alternative K_{e6} = <e6, e0>, eps w_1, eps w_2 are independent: WLOG
         w_1 = e'_1, w_2 = e'_3.
For general even n the same alternatives are used with x = e_{2k} (k < n/2) or x = e0
(k = n/2); run with --n 6 to validate against the n = 6 results.
Usage: python3 sym_p2.py TIME_LIMIT CASE [CASE ...] [--n 6]
"""
import json, os, subprocess, sys, tempfile, time
from itertools import combinations
import os, sys  # repository layout: shared modules are in ../../core
sys.path.insert(1, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "core"))
from sym_search import cases, normal_form, build_sym, fixed_dim
from sym_runner import write_cnf, KISSAT, DRAT
from du4perm_sat import extract, lut, check
n = int(sys.argv[sys.argv.index("--n") + 1]) if "--n" in sys.argv else 8
pairs = list(combinations(range(n), 2)); pidx = {pq: t for t, pq in enumerate(pairs)}
def zero(i, j):
    return [(1 << (pidx[(min(i, j), max(i, j))] * n + t), 0) for t in range(n)]
def desarg_v1():
    om = {0: [2], 2: [0, 2], 4: [6], 6: [4, 6]}
    rows = []
    idx = [0, 2, 4, 6]
    def vec(terms):
        out = []
        for t in range(n):
            m = 0
            for (i, j) in terms:
                if i != j:
                    m ^= 1 << (pidx[(min(i, j), max(i, j))] * n + t)
            out.append((m, 0))
        return out
    for a in idx:
        rows += vec([(a, j) for j in om[a]])
        for b in idx:
            if b > a:
                rows += vec([(a, j) for j in om[b]] + [(b, j) for j in om[a]])
    return rows
TMP = os.environ.get("SYM_TMP", tempfile.gettempdir())
log = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "logs", f"sym_p2_n{n}.jsonl")
TL = int(sys.argv[1])
args = [a for a in sys.argv[2:]]
if "--n" in args:
    i = args.index("--n"); del args[i:i + 2]
only_alt = None
if "--alt" in args:
    i = args.index("--alt"); only_alt = int(args[i + 1]); del args[i:i + 2]
args = [a for a in args if a != "--outnorm"]
for idx in map(int, args):
    p, la, A, lb, B = cases(n)[idx]
    ok, why, base = normal_form(p, A, B, n, la, lb)
    if not ok:
        print(idx, "excluded:", why, flush=True); continue
    k, kB = n - fixed_dim(A, n), n - fixed_dim(B, n)
    if 2 * k < n:
        x = 2 * k
        alts = [(f"K_e{x}=<e{x},e0>", zero(x, 0)), (f"K_e{x}=<e{x},e{x+1}>", zero(x, x + 1))]
        if "--outnorm" in sys.argv and n == 8 and k == 3 and kB == 3:
            def setv(i, j, target):
                return [(1 << (pidx[(min(i, j), max(i, j))] * n + t), target >> t & 1) for t in range(n)]
            alts = [("K_e6=<e6,e0>+out", zero(6, 0) + setv(6, 3, 1 << 1) + setv(6, 5, 1 << 3)),
                    ("K_e6=<e6,e7>+out", zero(6, 7) + setv(6, 1, 1 << 1) + setv(6, 3, 1 << 3) + setv(6, 5, 1 << 5))]
    elif n == 8 and kB == 2:
        alts = [("V1 closed, F4-normalised", desarg_v1())]
    else:
        alts = [("K_e0=<e0,e2>", zero(0, 2)), ("K_e0=<e0,e1>", zero(0, 1))]
    for ai, (name, rows) in enumerate(alts):
        if only_alt is not None and ai != only_alt:
            continue
        t0 = time.time()
        rec = {"case": idx, "A": la, "B": lb, "alt": name}
        try:
            enc, c, ell, dim = build_sym(n, A, B, True, False, base + rows)
        except ValueError as e:
            rec.update(status="UNSAT", note=str(e))
        else:
            cnf = os.path.join(TMP, f"p2_{idx}_{len(name)}_{os.getpid()}.cnf")
            write_cnf(enc, cnf, f"case {idx} {name}")
            r = subprocess.run([KISSAT, "-q", f"--time={TL}", cnf], capture_output=True, text=True)
            rec.update(bits=dim, vars=enc.pool.top, clauses=len(enc.clauses))
            if r.returncode == 10:
                model = [int(x) for line in r.stdout.splitlines() if line.startswith("v ")
                         for x in line[2:].split() if x != "0"]
                C, L = extract(model, c, ell, n)
                rec.update(status="SAT", check=check(lut(C, L, n), n), lut=lut(C, L, n))
            elif r.returncode == 20:
                rec["status"] = "UNSAT"
                if time.time() - t0 < 300:
                    proof = cnf + ".drat"
                    subprocess.run([KISSAT, "-q", "--no-binary", cnf, proof], capture_output=True)
                    v = subprocess.run([DRAT, cnf, proof], capture_output=True, text=True)
                    rec["drat"] = "VERIFIED" if "s VERIFIED" in v.stdout.replace("\r", "\n") else "NOT VERIFIED"
                    os.remove(proof)
            else:
                rec["status"] = "TIMEOUT"
            os.remove(cnf)
        rec.update(seconds=round(time.time() - t0, 1), when=time.strftime("%Y-%m-%d %H:%M:%S"))
        with open(log, "a") as f:
            f.write(json.dumps(rec) + "\n")
        print({a: b for a, b in rec.items() if a != "lut"}, flush=True)
