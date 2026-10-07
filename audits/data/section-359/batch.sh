#!/bin/bash
cd "$(dirname "$0")"
for v in A C; do for i in $(seq 0 47); do echo $v $i; done; done | xargs -P 4 -n 2 ./job.sh
echo BATCH-DONE $(ls out/*.txt | wc -l)
