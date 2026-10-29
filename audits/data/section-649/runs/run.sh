#!/bin/bash
R=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r649/fix; P=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r648
cd $R && xargs -P 4 -I{} bash -c "{}" < jobs.txt
( cd new && ./crem_interactions_validation > ../validation.txt 2>&1; echo EXIT $? >> ../validation.txt )
run(){ tag=$1; shift; mkdir -p reg/$tag; ( cd reg/$tag && $R/new/crem_interactions "$@" > out.txt 2> err.txt ); }
run 1 --mode statistical --phenomenon 1 --runs 16 --seed 466
run 2 --mode statistical --phenomenon 2 --runs 16 --seed 466
run L2 --mode statistical --phenomenon 2 --runs 16 --seed 466 --level 2
run 5Ps --mode statistical --phenomenon 5 --runs 50 --seed 551
run 5H --mode statistical --phenomenon 5 --runs 50 --seed 551 --pair proton,electron
run 6 --mode statistical --phenomenon 6 --runs 200 --seed 42
echo DONE
