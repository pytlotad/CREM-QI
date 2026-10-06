#!/bin/bash
cd "$(dirname "$0")"
{ for j in 0.25 0.75; do for v in A C; do for i in 0 1 2 3 4 5; do echo 2 $j $v $i; done; done; done
  for j in 0.1666666667 0.5 0.8333333333; do for i in 0 1 2 3; do echo 3 $j C $i; done; done; } | xargs -P 4 -n 4 ./job.sh
echo BATCH-DONE $(ls out/*.txt | wc -l)
