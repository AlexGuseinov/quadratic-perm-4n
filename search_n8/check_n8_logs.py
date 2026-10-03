"""Check the logs of the n = 8 search (run_n8.sh) for completeness.

Passes iff: both parts finished ("done in"), no "BUG" line, every level-3 index 0..N-1 appears
exactly once (as "done" or "pruned") with the hash of the level-3 list, the per-part COUNTS lines
show no node of dimension 6, 7 or 8, the solution files are empty, and the numbers of nodes
agree with expected_counts.json.
Usage: python3 search_n8/check_n8_logs.py [PREFIX]   (default <repository>/logs/crsub_n8_final)
       optional: --list list3_basis.txt  (compare the level-3 hashes with a separate listing)
"""
import json
import os
import re
import sys

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")   # repository root

prefix = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("--") else os.path.join(ROOT, "logs", "crsub_n8_final")
ok = True
seen = {}
totals = {}
for p in (0, 1):
    log = open(f"{prefix}_part{p}.log").read()
    if "BUG" in log:
        print(f"part {p}: BUG line present"); ok = False
    if "done in" not in log:
        print(f"part {p}: not finished"); ok = False
    for m in re.finditer(r"^L3 (\d+) W3 ([0-9a-f]{16}) (done|pruned)", log, re.M):
        i = int(m.group(1))
        if i in seen:
            print(f"index {i} appears twice"); ok = False
        if i % 2 != p:
            print(f"index {i} in the wrong part"); ok = False
        seen[i] = m.group(2)
    m = re.search(r"^COUNTS (.*)$", log, re.M)
    if not m:
        print(f"part {p}: no COUNTS line"); ok = False
    else:
        for kv in m.group(1).split():
            k, v = kv.split("=")
            totals[k] = totals.get(k, 0) + int(v)
    sol = open(f"{prefix}_sol_part{p}.txt").read().strip()
    if sol:
        print(f"part {p}: SOLUTIONS FOUND:\n{sol}"); ok = False
N = max(seen) + 1 if seen else 0
missing = [i for i in range(N) if i not in seen]
if missing:
    print("missing indices:", missing[:20]); ok = False
if "--list" in sys.argv:
    ref = {}
    for line in open(sys.argv[sys.argv.index("--list") + 1]):
        if line.startswith("L3"):
            f = line.split()
            ref[int(f[1])] = f[3]
    if ref.keys() != seen.keys() or any(ref[i] != seen[i] for i in ref):
        print("level-3 hashes differ from the listing"); ok = False
    else:
        print(f"level-3 hashes agree with the listing ({len(ref)} nodes)")
for d in (6, 7, 8):
    if totals.get(f"n{d}", 0):
        print(f"nodes of dimension {d}: {totals[f'n{d}']}"); ok = False
print(f"level-3 nodes: {N}; nodes by dimension (summed over parts; dimensions 1-2 are counted "
      f"in both parts): " + ", ".join(f"{d}: {totals.get(f'n{d}', 0)}" for d in range(1, 9)))
print("pruned (dim, cover, size) by dimension: " + "; ".join(
    f"{d}: {totals.get(f'pd{d}', 0)}, {totals.get(f'pc{d}', 0)}, {totals.get(f'ps{d}', 0)}" for d in range(3, 9)))
exp_file = os.path.join(ROOT, "expected_counts.json")
if os.path.exists(exp_file):
    exp = json.load(open(exp_file))["first_program"]
    for key, tag in (("nodes", "n"), ("returned_in_step_2", "ps"), ("returned_in_step_3", "pc")):
        for d, v in exp[key].items():
            if int(d) >= 3 and totals.get(f"{tag}{d}", 0) != v:
                print(f"{key}, dimension {d}: expected {v}, found {totals.get(f'{tag}{d}', 0)}")
                ok = False
    print("counts compared with expected_counts.json")
print("PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
