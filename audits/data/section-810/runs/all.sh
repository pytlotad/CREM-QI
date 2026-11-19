#!/bin/bash
# Audit 810: exact Sommerfeld remainder in the hydrogen level energy (default; CREM_NO_SOMMERFELD_EXACT=1 = before).  Base tree d478b489.
# old = git archive d478b489, new = tree, new_off = tree + CREM_NO_SOMMERFELD_EXACT=1.  Sequential, LC_ALL=C, memory cap (766).
set -u; export LC_ALL=C OPENBLAS_NUM_THREADS=1; ulimit -v 4000000
R=/home/teddy/Projekty/CREM-QI; S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s810r; D=$R/audits/data/section-810
rm -rf $S; mkdir -p $S/old $S/new
( cd $R && git archive d478b489 | tar -x -C $S/old ); ( cd $R && rsync -a --exclude .git --exclude audits ./ $S/new/ )
mkdir -p $S/new/audits/data/section-585; cp $R/audits/data/section-585/p585.cpp $S/new/audits/data/section-585/
F="-O2 -std=c++20 $(root-config --cflags) -I. -Imodules"; L="$(root-config --libs)"
for v in old new; do ( cd $S/$v && g++ $F audits/data/section-585/p585.cpp -o $S/p585_$v $L && g++ $F $D/p810.cpp -o $S/p810_$v $L ); [ -x $S/p585_$v ] && [ -x $S/p810_$v ] || { echo "BLAD kompilacja $v"; exit 1; }; done
echo "build probes $(date +%T)"
$S/p810_new > $D/out/p810_new.txt; CREM_NO_SOMMERFELD_EXACT=1 $S/p810_new > $D/out/p810_new_off.txt; $S/p810_old > $D/out/p810_old.txt
bash $R/audits/data/section-741/runs/probes.sh $S/p585_old $S/probes/old
bash $R/audits/data/section-741/runs/probes.sh $S/p585_new $S/probes/new
bash $R/audits/data/section-741/runs/probes.sh $S/p585_new $S/probes/new_off CREM_NO_SOMMERFELD_EXACT=1
for t in H2p H3p H4p H3s H3d Ps2p; do for w in new new_off; do s=SAME; cmp -s $S/probes/old/$t.out $S/probes/$w/$t.out && cmp -s $S/probes/old/$t.err $S/probes/$w/$t.err || s=DIFF; echo "probe $t old vs $w $s"; done; done > $D/out/probes810.txt
mkdir -p $D/out/probes_old $D/out/probes_new; cp $S/probes/old/*.err $D/out/probes_old/; cp $S/probes/new/*.err $D/out/probes_new/
echo "probes $(date +%T)"
( cd $S/old && make -j2 crem_interactions crem_interactions_validation >/dev/null 2>&1 ); ( cd $S/new && make -j2 crem_interactions crem_interactions_validation >/dev/null 2>&1 )
for v in old new; do [ -x $S/$v/crem_interactions ] && [ -x $S/$v/crem_interactions_validation ] || { echo "BLAD make $v"; exit 1; }; done
echo "build $(date +%T)"
for v in old new; do ( cd $S/$v && ./crem_interactions_validation > $D/out/validation_$v.txt 2>&1; echo EXIT $? >> $D/out/validation_$v.txt ); echo "validation $v $(date +%T)"; done
run(){ v=$1; tag=$2; shift 2; mkdir -p $S/reg/$v/$tag; ( cd $S/reg/$v/$tag && $S/${v%%_*}/crem_interactions "$@" > out.txt 2> err.txt ); }
for v in old new new_off; do
  if [ $v = new_off ]; then export CREM_NO_SOMMERFELD_EXACT=1; else unset CREM_NO_SOMMERFELD_EXACT; fi
  run $v 1 --mode statistical --phenomenon 1 --runs 16 --seed 466
  run $v 2 --mode statistical --phenomenon 2 --runs 16 --seed 466
  run $v L2 --mode statistical --phenomenon 2 --runs 16 --seed 466 --level 2
  run $v 5Ps --mode statistical --phenomenon 5 --runs 50 --seed 551
  run $v 5H --mode statistical --phenomenon 5 --runs 50 --seed 551 --pair proton,electron
  run $v 6 --mode statistical --phenomenon 6 --runs 200 --seed 42
  echo "reg $v $(date +%T)"
done
unset CREM_NO_SOMMERFELD_EXACT
for w in new new_off; do for tag in 1 2 L2 5Ps 5H 6; do s=SAME; for x in out.txt err.txt; do cmp -s $S/reg/old/$tag/$x $S/reg/$w/$tag/$x || s=DIFF; done; echo "old vs $w $tag $s"; done; done > $D/out/regression810.txt
mkdir -p $D/out/reg5H; cp $S/reg/old/5H/err.txt $D/out/reg5H/old.err; cp $S/reg/new/5H/err.txt $D/out/reg5H/new.err
echo "ALLDONE $(date +%T)"
