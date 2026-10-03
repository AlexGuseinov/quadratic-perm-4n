#!/bin/bash
# Solve all 2^DEPTH cubes, P at a time; each cube CNF is created on the fly and deleted.
# Usage (from earlier_approaches/sat): bash vps/run_cubes.sh DEPTH [P]
DEPTH=${1:-10}; P=${2:-$(nproc)}
mkdir -p results tmpcnf
seq 0 $(( (1 << DEPTH) - 1 )) | xargs -P "$P" -I{} bash -c '
  m={}; out=results/cube_${m}.out
  if [ -f $out ] && grep -q "^s " $out; then exit 0; fi
  python3 vps/cube_cnf.py '"$DEPTH"' $m tmpcnf/c$m.cnf
  ~/tools/kissat/build/kissat -q tmpcnf/c$m.cnf > $out 2>&1
  rm -f tmpcnf/c$m.cnf'
grep -h "^s " results/*.out | sort | uniq -c
