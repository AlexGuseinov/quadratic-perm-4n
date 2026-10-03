#!/bin/bash
# one-line status of the n = 8 runs
cd "$(dirname "$0")/../.."
for f in logs/kissat_d8_sb3.log logs/kissat_d8_noperm.log; do
  s=$(grep '^s ' $f)
  t=$(grep -E '^c [-a-z.] +[0-9]' $f | tail -1 | awk '{print $3"s"}')
  echo "$(date +%H:%M) $f ${s:-running} $t"
done | tee -a logs/poll.log
ps -eo pid,etime,args | grep -cE "^ *[0-9]+ +[0-9:-]+ bin/kissat"
