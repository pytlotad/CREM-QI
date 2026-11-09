#!/bin/bash
# Audit 768: rho3 table k = 2 in CREM_LS_FROM_PRECESSION (opt-in).  Base tree d80690a5 (767) vs working tree.  Sequential.  LC_ALL=C.  Memory cap (766).
set -u; export LC_ALL=C OPENBLAS_NUM_THREADS=1; ulimit -v 4000000
R=/home/teddy/Projekty/CREM-QI; S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s768r; D=$R/audits/data/section-768
rm -rf $S; mkdir -p $S/old $S/new
( cd $R && git archive d80690a5 | tar -x -C $S/old ); ( cd $R && rsync -a --exclude .git --exclude audits ./ $S/new/ )
mkdir -p $S/new/audits/data/section-585; cp $R/audits/data/section-585/p585.cpp $S/new/audits/data/section-585/
F="-O2 -std=c++20 $(root-config --cflags) -I. -Imodules"; L="$(root-config --libs)"
for v in old new; do ( cd $S/$v && g++ $F audits/data/section-585/p585.cpp -o $S/p585_$v $L && g++ $F $R/audits/data/section-763/p763.cpp -o $S/p763_$v $L ); [ -x $S/p585_$v ] && [ -x $S/p763_$v ] || { echo "BLAD kompilacja $v"; exit 1; }; done
$S/p763_new > $D/out/p763_default.txt; CREM_LS_FROM_PRECESSION=1 $S/p763_new > $D/out/p763_switch.txt; CREM_LS_FROM_PRECESSION=1 $S/p763_old > $D/out/p763_old_switch.txt
for v in old new; do bash $R/audits/data/section-741/runs/probes.sh $S/p585_$v $S/probes/$v; done
bash $R/audits/data/section-741/runs/probes.sh $S/p585_new $S/probes/sw CREM_LS_FROM_PRECESSION=1
for t in H2p H3p H4p H3s H3d Ps2p; do s=SAME; cmp -s $S/probes/old/$t.out $S/probes/new/$t.out && cmp -s $S/probes/old/$t.err $S/probes/new/$t.err || s=DIFF
  w=SAME; cmp -s $S/probes/new/$t.err $S/probes/sw/$t.err || w=DIFF; echo "probe $t old vs new $s; new vs switch $w"; done > $D/out/probes768.txt
mkdir -p $D/out/probes_switch $D/out/probes_default; cp $S/probes/sw/*.err $D/out/probes_switch/; cp $S/probes/new/*.err $D/out/probes_default/
( cd $S/old && make -j2 crem_interactions >/dev/null 2>&1 ); ( cd $S/new && make -j2 crem_interactions crem_interactions_validation >/dev/null 2>&1 )
[ -x $S/old/crem_interactions ] && [ -x $S/new/crem_interactions ] && [ -x $S/new/crem_interactions_validation ] || { echo "BLAD make"; exit 1; }
( cd $S/new && ./crem_interactions_validation > $D/out/validation.txt 2>&1; echo EXIT $? >> $D/out/validation.txt )
run(){ v=$1; tag=$2; shift 2; mkdir -p $S/reg/$v/$tag; ( cd $S/reg/$v/$tag && $S/$v/crem_interactions "$@" > out.txt 2> err.txt ); }
for v in old new; do
  run $v 1 --mode statistical --phenomenon 1 --runs 16 --seed 466
  run $v 2 --mode statistical --phenomenon 2 --runs 16 --seed 466
  run $v L2 --mode statistical --phenomenon 2 --runs 16 --seed 466 --level 2
  run $v 5Ps --mode statistical --phenomenon 5 --runs 50 --seed 551
  run $v 5H --mode statistical --phenomenon 5 --runs 50 --seed 551 --pair proton,electron
  run $v 6 --mode statistical --phenomenon 6 --runs 200 --seed 42
done
for tag in 1 2 L2 5Ps 5H 6; do s=SAME; for x in out.txt err.txt; do cmp -s $S/reg/old/$tag/$x $S/reg/new/$tag/$x || s=DIFF; done; echo "old vs new $tag $s"; done > $D/out/regression768.txt
echo ALLDONE
