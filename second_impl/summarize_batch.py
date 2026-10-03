"""Summary of the batch run of the second implementation over all roots (batch_results_*.txt).

Checks that every root index of the listing appears exactly once, that no root has solutions,
and prints the totals of nodes and pruned nodes by dimension (3..8) and the total time.
The totals are compared with ../expected_counts.json (if present).
Usage: python3 summarize_batch.py [RESULTS ...]   (default: batch_results_0.txt batch_results_1.txt)
"""
import json
import os
import re
import sys

here = os.path.dirname(os.path.abspath(__file__))
files = sys.argv[1:] or [os.path.join(here, f"batch_results_{i}.txt") for i in (0, 1)]
listing = os.path.join(os.environ.get("N8_DATA", os.path.join(here, "data")), "list3_basis.txt")
roots = {int(l.split()[1]) for l in open(listing) if l.startswith("L3")}
seen, tot, secs, ok = {}, {}, 0.0, True
for f in files:
    for line in open(f):
        if not line.startswith("ROOT"):
            continue
        kv = dict(re.findall(r"(\w+)=(\S+)", line))
        i = int(kv["idx"])
        if i in seen:
            print("duplicate root", i); ok = False
        seen[i] = kv
        if kv["sols"] != "0":
            print("SOLUTIONS at root", i, line.strip()); ok = False
        for key in ("nodes", "pdim", "psize", "pcov"):
            vals = [int(x) for x in kv[key].split(",")]
            t = tot.setdefault(key, [0] * len(vals))
            for k, v in enumerate(vals):
                t[k] += v
        secs += float(kv["secs"])
missing = sorted(roots - seen.keys())
if missing:
    print("missing roots:", missing[:20], "..." if len(missing) > 20 else ""); ok = False
print(f"roots: {len(seen)} of {len(roots)}; total time {secs / 3600:.2f} h")
for key in ("nodes", "pdim", "psize", "pcov"):
    print(f"{key:6s} by dimension 3..8: {tot.get(key)}")
exp_file = os.path.join(here, "..", "expected_counts.json")
if os.path.exists(exp_file):
    exp = json.load(open(exp_file))["second_program"]
    got = {"nodes": tot.get("nodes"), "returned_in_step_2": tot.get("psize"),
           "returned_in_step_3": tot.get("pcov")}
    for key, arr in got.items():
        for d, v in exp[key].items():
            if arr is None or arr[int(d) - 3] != v:
                print(f"{key}, dimension {d}: expected {v}, found {None if arr is None else arr[int(d) - 3]}")
                ok = False
    if tot.get("pdim") and any(tot["pdim"]):
        print("nodes with dim P_a > 2 occurred"); ok = False
    print("counts compared with expected_counts.json")
print("PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
