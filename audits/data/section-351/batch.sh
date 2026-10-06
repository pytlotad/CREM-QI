#!/bin/bash
cd "$(dirname "$0")"; mkdir -p out
echo 3 15 22 24 38 36 34 0 | tr ' ' '\n' | xargs -P 4 -n 1 ./job.sh
echo BATCH-DONE $(ls out/*.txt | wc -l)
