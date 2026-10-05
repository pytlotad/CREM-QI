#!/bin/bash
cd "$(dirname "$0")"; mkdir -p out
for v in A C; do for i in $(seq 0 47); do echo $v $i; done; done | xargs -P 4 -n 2 ./job.sh
echo BATCH-DONE $(ls out/*.txt | wc -l)
cd /home/teddy/Projekty/Positronium && make validation-small > "$OLDPWD/walidacja-small.txt" 2>&1; git checkout -- distributions/ 2>/dev/null
tail -2 "$OLDPWD/walidacja-small.txt"
