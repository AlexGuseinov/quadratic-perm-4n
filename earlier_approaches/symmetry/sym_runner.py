"""Run every symmetry case of sym_search.py with kissat.

For each case: write the CNF, run kissat with a time limit; if SAT, rebuild the
encoding, read the model and check the LUT independently (bijective, DU = 4);
if UNSAT and the run was short, re-run with a DRAT proof and check it with
drat-trim.  One JSON line per case in logs/sym_n{n}.jsonl.
Usage: python3 sym_runner.py n TIME_LIMIT [--workers 2] [--only 3,5,7] [--primes 3,5]
                            [--proof-max 120] [--noperm] [--reduce]
  --noperm: drop the permutation condition (is there any symmetric beta at all?)
  --reduce: keep rank (and, with --noperm, witness) constraints for orbit representatives only
  --nf: add the lemma checks and normal-form rows of sym_search.normal_form (N1-N10, SB1', SB2', SB2'', SB3')
"""
import json
import os
import subprocess
import sys
import tempfile
import time
from multiprocessing import Pool

import os, sys  # repository layout: shared modules are in ../../core
sys.path.insert(1, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "core"))
from sym_search import cases, build_sym, equivariant_basis, normal_form
from du4perm_sat import extract, lut, check

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")   # repository root (logs/ is there)
KISSAT = os.path.join(HERE, "..", "bin", "kissat")
DRAT = os.path.join(HERE, "..", "bin", "drat-trim")
TMP = os.environ.get("SYM_TMP", tempfile.gettempdir())


def write_cnf(enc, path, head):
    with open(path, "w") as f:
        f.write(f"c {head}\n")
        f.write(f"p cnf {enc.pool.top} {len(enc.clauses)}\n")
        for cl in enc.clauses:
            f.write(" ".join(map(str, cl)) + " 0\n")


def run_case(args):
    n, idx, tlimit, proof_max, perm, reduce, nf = args
    p, la, A, lb, B = cases(n)[idx]
    rows = None
    if nf:
        ok, why, rows = normal_form(p, A, B, n, la, lb)
        if not ok:
            return {"n": n, "case": idx, "p": p, "A": la, "B": lb, "status": "UNSAT",
                    "seconds": 0.0, "note": why, "when": time.strftime("%Y-%m-%d %H:%M:%S")}
    try:
        enc, c, ell, dim = build_sym(n, A, B, perm, reduce, rows)
    except ValueError as e:
        # a constraint became an empty clause during encoding (e.g. beta = 0)
        dim = len(equivariant_basis(A, B, n)[0])
        return {"n": n, "case": idx, "p": p, "A": la, "B": lb, "free_beta_bits": dim,
                "status": "UNSAT", "seconds": 0.0,
                "note": "normal-form equations inconsistent" if "normal-form" in str(e)
                        else "empty clause during encoding",
                "when": time.strftime("%Y-%m-%d %H:%M:%S")}
    head = f"n={n} case={idx} p={p} A={la} B={lb} free_beta_bits={dim}"
    cnf = os.path.join(TMP, f"sym_{n}_{idx}_{int(perm)}{int(reduce)}{int(nf)}_{os.getpid()}.cnf")
    write_cnf(enc, cnf, head)
    t0 = time.time()
    r = subprocess.run([KISSAT, "-q", f"--time={tlimit}", cnf],
                       capture_output=True, text=True)
    dt = time.time() - t0
    rec = {"n": n, "case": idx, "p": p, "A": la, "B": lb, "free_beta_bits": dim,
           "perm": perm, "reduce": reduce, "nf": nf, "vars": enc.pool.top, "clauses": len(enc.clauses),
           "seconds": round(dt, 1), "time_limit": tlimit,
           "kissat": "4.0.4", "when": time.strftime("%Y-%m-%d %H:%M:%S")}
    if r.returncode == 10:
        model = []
        for line in r.stdout.splitlines():
            if line.startswith("v "):
                model += [int(x) for x in line[2:].split() if x != "0"]
        C, L = extract(model, c, ell, n)
        table = lut(C, L, n)
        ok, du = check(table, n)
        rec.update(status="SAT", check_perm=ok, check_du=du, lut=table,
                   beta={f"{i},{j}": v for (i, j), v in C.items()})
    elif r.returncode == 20:
        rec["status"] = "UNSAT"
        if dt <= proof_max:
            proof = cnf + ".drat"
            subprocess.run([KISSAT, "-q", "--no-binary", cnf, proof],
                           capture_output=True, text=True)
            v = subprocess.run([DRAT, cnf, proof], capture_output=True, text=True)
            out = v.stdout.replace("\r", "\n")
            rec["drat"] = "VERIFIED" if "s VERIFIED" in out else "NOT VERIFIED"
            try:
                os.remove(proof)
            except OSError:
                pass
    else:
        rec["status"] = "TIMEOUT"
    os.remove(cnf)
    return rec


def main():
    n, tlimit = int(sys.argv[1]), int(sys.argv[2])
    workers = int(sys.argv[sys.argv.index("--workers") + 1]) if "--workers" in sys.argv else 2
    proof_max = int(sys.argv[sys.argv.index("--proof-max") + 1]) if "--proof-max" in sys.argv else 120
    cs = cases(n)
    idxs = list(range(len(cs)))
    if "--only" in sys.argv:
        idxs = [int(x) for x in sys.argv[sys.argv.index("--only") + 1].split(",")]
    if "--primes" in sys.argv:
        ps = {int(x) for x in sys.argv[sys.argv.index("--primes") + 1].split(",")}
        idxs = [i for i in idxs if cs[i][0] in ps]
    # easiest first: fewest free beta bits
    idxs.sort(key=lambda i: len(equivariant_basis(cs[i][2], cs[i][4], n)[0]))
    if "--order" in sys.argv:
        idxs = [int(x) for x in sys.argv[sys.argv.index("--order") + 1].split(",")]
    perm = "--noperm" not in sys.argv
    reduce = "--reduce" in sys.argv
    nf = "--nf" in sys.argv
    tag = ("" if perm else "_noperm") + ("_red" if reduce else "") + ("_nf" if nf else "")
    out = os.path.join(ROOT, "logs", f"sym_n{n}{tag}.jsonl")
    with Pool(workers) as pool:
        for rec in pool.imap_unordered(run_case, [(n, i, tlimit, proof_max, perm, reduce, nf) for i in idxs]):
            with open(out, "a") as f:
                f.write(json.dumps(rec) + "\n")
            extra = f" [{rec['note']}]" if "note" in rec else ""
            if rec["status"] == "SAT":
                extra = f" perm={rec['check_perm']} DU={rec['check_du']}"
            if "drat" in rec:
                extra += f" drat={rec['drat']}"
            print(f"case {rec['case']:3d} p={rec['p']:<3d} A={rec['A']:<12s} B={rec['B']:<12s} "
                  f"bits={rec.get('free_beta_bits', -1):3d} {rec['status']:7s} {rec['seconds']:8.1f}s{extra}",
                  flush=True)


if __name__ == "__main__":
    main()
