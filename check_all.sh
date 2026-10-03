#!/bin/bash
# Quick verification of all main claims (about 1 minute on a laptop).
# Run from the repository folder:  bash check_all.sh
set -e
cd "$(dirname "$0")"

echo "== 1. Basics: rank counts, x^5 on GF(64), n=4 constant-rank bound"
echo "   expect: match=True; DU: 4; component ranks {4: 63}; maximal dimension: 3"
python3 core/validate_basics.py

echo; echo "== 2. MacWilliams test (L5)"
echo "   expect: n= 4 ... d= 4: IMPOSSIBLE ;  n= 8 ... d= 8: consistent"
python3 core/macwilliams_alt.py | grep -E "match|n= 4 \(n mod 4 = 0\) d= 4|n= 8 \(n mod 4 = 0\) d= 8"

echo; echo "== 3. Exact SAT model D(n)"
echo "   expect: n=4 UNSAT;  n=6 SAT with permutation=True DU=4"
python3 core/du4perm_sat.py 4 --sb3 | tail -1
python3 core/du4perm_sat.py 6 --sb3 | sed -n 2p

echo; echo "== 4. n=8 with Desarguesian kernel spread (Theorem L7)"
echo "   expect: UNSAT (also without permutation condition)"
python3 core/du4perm_sat.py 8 --sb3 --desarg | tail -1
python3 core/du4perm_sat.py 8 --desarg --noperm | tail -1

echo; echo "== 5. Hand-checkable examples for L3 and L7"
echo "   expect: [1] True ...  [2] permutation=True, DU=4, omegas found=2, Phi F_4-bilinear=True"
python3 core/check_x5.py

echo; echo "All quick checks done."
