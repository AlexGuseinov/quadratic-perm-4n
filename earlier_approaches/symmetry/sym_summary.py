"""Summary of the symmetric-case computation for n = 8 (Theorem sym8): for each of the
93 pairs (A, B) print how it was excluded, from the logs written by sym_runner.py
(--nf), sym_twostage.py (--nf), sym_case44.py, sym_case32.py and sym_p2.py."""
import hashlib, json, os
from sym_search import cases, normal_form, affine_space
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")   # repository root (logs/ is there)
def load(name):
    path = os.path.join(ROOT, "logs", name)
    return [json.loads(l) for l in open(path)] if os.path.exists(path) else []
n = 8
files = ["sym_search.py", "sym_runner.py", "sym_twostage.py", "sym_case32.py", "sym_case44.py",
         "sym_p2.py", "beta_enum.c", "../../core/du4perm_sat.py"]
for fn in files:
    h = hashlib.sha256(open(os.path.join(HERE, fn), "rb").read()).hexdigest()[:16]
    print(f"# {os.path.basename(fn):18s} sha256 {h}")
runner = {}
for r in load("sym_n8_nf.jsonl"):
    if r["status"] == "UNSAT":
        runner[r["case"]] = r
enum = {r["case"]: r for r in load("sym2_n8_nf.jsonl") if r["status"] == "UNSAT"}
c44 = load("sym_case44.jsonl"); c32 = load("sym_case32.jsonl"); p2 = load("sym_p2_n8.jsonl")
def all_unsat(recs):
    return bool(recs) and all(r["status"] == "UNSAT" for r in recs)
open_cases = []
for i, (p, la, A, lb, B) in enumerate(cases(n)):
    ok, why, rows = normal_form(p, A, B, n, la, lb)
    if not ok:
        how = "Lemma " + why
    elif i == 44 and all_unsat(c44) and {r["v"] for r in c44} == {2, 4, 6, 8, 10, 12, 14, 16}:
        how = "8 alternatives, exact enumeration / inconsistent (sym_case44.py)"
    elif (i == 32 and all_unsat(c32) and {r["class"] for r in c32} == set(range(42))
          and all(r.get("drat") == "VERIFIED" for r in c32)):
        how = "42 classes, kissat UNSAT + DRAT (sym_case32.py)"
    elif i in (11, 12, 16):
        expected = {11: {"K_e6=<e6,e0>", "K_e6=<e6,e7>"},
                    12: {"K_e6=<e6,e0>+out", "K_e6=<e6,e7>+out"},
                    16: {"V1 closed, F4-normalised"}}[i]
        rr = [r for r in p2 if r["case"] == i and r["alt"] in expected
              and r["status"] == "UNSAT" and r.get("drat") == "VERIFIED"]
        if {r["alt"] for r in rr} == expected:
            how = f"{len(expected)} alternative(s), kissat UNSAT + DRAT (sym_p2.py)"
        else:
            how = "OPEN"
            open_cases.append(i)
    elif i in runner and runner[i].get("drat") == "VERIFIED":
        how = f"kissat UNSAT + DRAT ({runner[i]['seconds']} s)"
    elif i in enum:
        how = f"exact enumeration: {enum[i].get('stage1_betas', 0)} admissible beta"
    elif i in runner:
        how = "UNSAT (" + runner[i].get("note", "kissat") + ")"
    else:
        how = "OPEN"
        open_cases.append(i)
    print(f"{i:3d}  p={p:<3d} A={la:<12s} B={lb:<12s} {how}")
print("open:", open_cases if open_cases else "none")
