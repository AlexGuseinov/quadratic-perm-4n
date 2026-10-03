#!/bin/sh
# Theorem thm:n8 (no quadratic DU-4 permutation of F_2^8): the complete search, in two parallel
# parts, with a provenance record (compiler, command lines, sha256 of sources, binary, rank table).
# Usage: sh run_n8.sh [RANK6_TABLE]     (default good8.bin; built by the program if missing)
# Afterwards: python3 check_n8_logs.py  (checks that all 486 level-3 nodes were processed)
set -e
cd "$(dirname "$0")"
CFLAGS="-O3 -march=native"
gcc $CFLAGS -o crsub_enum crsub_enum.c -lm
GOOD=${1:-good8.bin}
[ -f "$GOOD" ] || ./crsub_enum 8 --good "$GOOD" --list3 > /dev/null 2>&1
P=logs/crsub_n8_final_provenance.txt
{
  echo "started: $(date -Iseconds)"
  uname -a
  gcc --version | head -1
  echo "build: gcc $CFLAGS -o crsub_enum crsub_enum.c -lm"
  sha256sum crsub_enum.c crsub_common.h crsub_enum "$GOOD"
} > $P
for p in 0 1; do
  echo "run: ./crsub_enum 8 --good $GOOD --part $p 2 --out logs/crsub_n8_final_sol_part$p.txt" >> $P
  ./crsub_enum 8 --good "$GOOD" --part $p 2 --out logs/crsub_n8_final_sol_part$p.txt \
      2> logs/crsub_n8_final_part$p.log &
done
wait
echo "finished: $(date -Iseconds)" >> $P
