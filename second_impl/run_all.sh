#!/bin/sh
# Verification of Theorem thm:n8 that does not use the first program (crsub_enum.c).
#
# Inputs taken from the repository (both produced by the first program, both CHECKED here in
# step 2 against the definitions, so they need not be trusted):
#   data/cert8.txt.xz     certificate of the symmetry reduction at dimensions 1 and 2
#   data/list3_basis.txt  the 486 nodes of dimension 3
# Steps:
#   1. task1_pfaffian.py   table of the forms of rank 6 by Pfaffians (count checked against the
#                          closed formula); writes data/rank6_pf.bin             (seconds)
#   2. verify_cert.py      checks the certificate and the list of the 486 nodes  (about 1 minute)
#   3. search.py --batch   searches the 486 subtrees without symmetry reduction
#                          (about 13 hours of process time in total, 0.6 GB per process)
#   4. summarize_batch.py  all 486 nodes processed, no solution, no node of dimension 6, node
#                          counts equal to ../expected_counts.json; prints PASS or FAIL
# Usage: sh run_all.sh [PROCESSES]      (default 2; step 3 can be interrupted and resumed)
# Needs Python 3.11 or later with NumPy 2.
set -e
cd "$(dirname "$0")"
P=${1:-2}
python3 task1_pfaffian.py
python3 verify_cert.py
i=0
while [ "$i" -lt "$P" ]; do
  python3 search.py --batch data/list3_basis.txt --out rerun_results_$i.txt --start "$i" --step "$P" \
      > rerun_batch$i.log 2>&1 &
  i=$((i + 1))
done
wait
python3 summarize_batch.py rerun_results_*.txt
