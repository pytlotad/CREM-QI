#!/bin/bash
# Audit 826: photon target on the whole rung (mechanical path; CREM_ACTION_RUNG_RELATIVE=1 = before).  Base tree cbd40b65 (825).
# (1) cascades: set of 825 on the NEW tree (48 runs, P3); (2) new_off on the 12 4f runs (byte check against 825); (3) probes 741
# old / new / new_off; (4) validation old / new; (5) regression old / new / new_off (sequential).  LC_ALL=C, memory caps (766).
set -u; export LC_ALL=C OPENBLAS_NUM_THREADS=1; ulimit -v 4000000
R=/home/teddy/Projekty/CREM-QI; S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s826; D=$R/audits/data/section-826
B0=cbd40b65; SW=CREM_ACTION_RUNG_RELATIVE=1
rm -rf $S; mkdir -p $S/old $S/new $D/out
( cd $R && git archive $B0 -- . ':(exclude)audits' | tar -x -C $S/old ); ( cd $R && rsync -a --exclude .git --exclude audits ./ $S/new/ )
for v in old new; do mkdir -p $S/$v/audits/data/section-585 $S/$v/audits/data/section-787
  git -C $R show $B0:audits/data/section-585/p585.cpp > $S/$v/audits/data/section-585/p585.cpp
  git -C $R show $B0:audits/data/section-787/p787.cpp > $S/$v/audits/data/section-787/p787.cpp; done
F="-O2 -std=c++20 $(root-config --cflags) -I. -Imodules"; L="$(root-config --libs)"
for v in old new; do ( cd $S/$v && g++ $F audits/data/section-585/p585.cpp -o $S/p585_$v $L ); [ -x $S/p585_$v ] || { echo "BLAD kompilacja p585 $v"; exit 1; }; done
# p787 exactly as 825 built it (-O3 -march=native), from the new tree
( cd $S/new && g++ -O3 -march=native -std=c++20 $(root-config --cflags) -I. -Imodules audits/data/section-787/p787.cpp -o $S/p787 $(root-config --libs) ) || { echo "BLAD kompilacja p787"; exit 1; }
echo "build probes $(date +%T)"
# (1) + (2) cascades
O=$D/out/runs; mkdir -p $O/new $O/off; export S O SW
{ for sys in "H0 PAIR=proton,electron PHEN=1" "H1 PAIR=proton,electron PHEN=2" "pPs PHEN=1" "oPs PHEN=2"; do set -- $sys; tag=$1; shift
    for st in "4f 4 3 30" "4d 4 2 30" "4p 4 1 30" "3d 3 2 40"; do set -- $st $*; s=$1 n=$2 l=$3 obs=$4; shift 4
      for i in 0 1 2; do echo "new/${tag}_${s}_$i $i NLEV=$n LHBAR=$l OBS=$obs $*"
        [ $s = 4f ] && echo "off/${tag}_${s}_$i $i NLEV=$n LHBAR=$l OBS=$obs $* $SW"; done; done; done; } |
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; seed=$1; shift
  env CREM_MECHANICAL_HAZARD_SCALE=8e5 CREM_DEBUG_LADDER=1 CREM_DEBUG_LADDER_SPLIT=1 CREM_PHASE_TRACE=1 "$@" $S/p787 $seed >$O/$name.out 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"'
echo "cascades $(date +%T)"
# (3) probes
bash $R/audits/data/section-741/runs/probes.sh $S/p585_old $S/probes/old
bash $R/audits/data/section-741/runs/probes.sh $S/p585_new $S/probes/new
bash $R/audits/data/section-741/runs/probes.sh $S/p585_new $S/probes/new_off $SW
for t in H2p H3p H4p H3s H3d Ps2p; do for w in new new_off; do s=SAME; cmp -s $S/probes/old/$t.out $S/probes/$w/$t.out && cmp -s $S/probes/old/$t.err $S/probes/$w/$t.err || s=DIFF; echo "probe $t old vs $w $s"; done; done > $D/out/probes826.txt
mkdir -p $D/out/probes_old $D/out/probes_new; cp $S/probes/old/*.err $D/out/probes_old/; cp $S/probes/new/*.err $D/out/probes_new/
echo "probes $(date +%T)"
# (4) validation
( cd $S/old && make -j2 crem_interactions crem_interactions_validation >/dev/null 2>&1 ); ( cd $S/new && make -j2 crem_interactions crem_interactions_validation >/dev/null 2>&1 )
for v in old new; do [ -x $S/$v/crem_interactions ] && [ -x $S/$v/crem_interactions_validation ] || { echo "BLAD make $v"; exit 1; }; done
for v in old new; do ( cd $S/$v && ./crem_interactions_validation > $D/out/validation_$v.txt 2>&1; echo EXIT $? >> $D/out/validation_$v.txt ); echo "validation $v $(date +%T)"; done
# (5) regression
run(){ v=$1; tag=$2; shift 2; mkdir -p $S/reg/$v/$tag; t0=$(date +%s); ( cd $S/reg/$v/$tag && $S/${v%%_*}/crem_interactions "$@" > out.txt 2> err.txt ); echo "$v $tag $(( $(date +%s) - t0 ))" >> $D/out/walltime826.txt; }
: > $D/out/walltime826.txt
for v in old new new_off; do
  if [ $v = new_off ]; then export CREM_ACTION_RUNG_RELATIVE=1; else unset CREM_ACTION_RUNG_RELATIVE; fi
  run $v 1 --mode statistical --phenomenon 1 --runs 16 --seed 826
  run $v 2 --mode statistical --phenomenon 2 --runs 16 --seed 826
  run $v L2 --mode statistical --phenomenon 2 --runs 16 --seed 826 --level 2
  run $v 5Ps --mode statistical --phenomenon 5 --runs 50 --seed 8261
  run $v 5H --mode statistical --phenomenon 5 --runs 50 --seed 8261 --pair proton,electron
  run $v 6 --mode statistical --phenomenon 6 --runs 200 --seed 8262
  echo "reg $v $(date +%T)"
done
unset CREM_ACTION_RUNG_RELATIVE
for w in new new_off; do for tag in 1 2 L2 5Ps 5H 6; do s=SAME; for x in out.txt err.txt; do cmp -s $S/reg/old/$tag/$x $S/reg/$w/$tag/$x || s=DIFF; done; echo "old vs $w $tag $s"; done; done > $D/out/regression826.txt
mkdir -p $D/out/reg; for v in old new; do for tag in 1 2 L2 5Ps 5H 6; do cp $S/reg/$v/$tag/out.txt $D/out/reg/${v}_$tag.out; done; done
echo "ALLDONE $(date +%T)"
