#!/bin/bash
R=$R658
( cd $R/new && ./crem_interactions_validation > ../validation.txt 2>&1; echo EXIT $? >> ../validation.txt )
run(){ v=$1; tag=$2; shift 2; b=$v; e=""; [ $v = flag ] && { b=new; e="CREM_H_HAMILTONIAN_TRANSPORT=1"; }; mkdir -p $R/reg/$v/$tag; ( cd $R/reg/$v/$tag && env $e $R/$b/crem_interactions "$@" > out.txt 2> err.txt ); }
for v in old new flag; do
run $v 1 --mode statistical --phenomenon 1 --runs 16 --seed 466
run $v 2 --mode statistical --phenomenon 2 --runs 16 --seed 466
run $v L2 --mode statistical --phenomenon 2 --runs 16 --seed 466 --level 2
run $v 5Ps --mode statistical --phenomenon 5 --runs 50 --seed 551
run $v 5H --mode statistical --phenomenon 5 --runs 50 --seed 551 --pair proton,electron
run $v 6 --mode statistical --phenomenon 6 --runs 200 --seed 42
done
echo DONE
