#!/bin/bash
# Audit 818: default Plummer floor scale 0.05 -> 0.01 r* (CREM_SEPARATION_FLOOR_SCALE=0.05 = before).  Base tree c2e76594.
# old = git archive c2e76594, new = tree, new_off = tree + CREM_SEPARATION_FLOOR_SCALE=0.05.  Sequential, LC_ALL=C, memory cap (766).
set -u; export LC_ALL=C OPENBLAS_NUM_THREADS=1; ulimit -v 4000000
R=/home/teddy/Projekty/CREM-QI; S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s818; D=$R/audits/data/section-818
rm -rf $S; mkdir -p $S/old $S/new $D/out
( cd $R && git archive c2e76594 -- . ':(exclude)audits' | tar -x -C $S/old ); ( cd $R && rsync -a --exclude .git --exclude audits ./ $S/new/ )
for v in old new; do mkdir -p $S/$v/audits/data/section-585; git -C $R show c2e76594:audits/data/section-585/p585.cpp > $S/$v/audits/data/section-585/p585.cpp; done
F="-O2 -std=c++20 $(root-config --cflags) -I. -Imodules"; L="$(root-config --libs)"
for v in old new; do ( cd $S/$v && g++ $F audits/data/section-585/p585.cpp -o $S/p585_$v $L ); [ -x $S/p585_$v ] || { echo "BLAD kompilacja $v"; exit 1; }; done
( cd $S/new && g++ $F $R/audits/data/section-814/p814.cpp -o $S/p814_new $L ); [ -x $S/p814_new ] || { echo "BLAD kompilacja p814"; exit 1; }
echo "build probes $(date +%T)"
bash $R/audits/data/section-741/runs/probes.sh $S/p585_old $S/probes/old
bash $R/audits/data/section-741/runs/probes.sh $S/p585_new $S/probes/new
bash $R/audits/data/section-741/runs/probes.sh $S/p585_new $S/probes/new_off CREM_SEPARATION_FLOOR_SCALE=0.05
for t in H2p H3p H4p H3s H3d Ps2p; do for w in new new_off; do s=SAME; cmp -s $S/probes/old/$t.out $S/probes/$w/$t.out && cmp -s $S/probes/old/$t.err $S/probes/$w/$t.err || s=DIFF; echo "probe $t old vs $w $s"; done; done > $D/out/probes818.txt
mkdir -p $D/out/probes_old $D/out/probes_new; cp $S/probes/old/*.err $D/out/probes_old/; cp $S/probes/new/*.err $D/out/probes_new/
echo "probes $(date +%T)"
( cd $S/old && make -j2 crem_interactions crem_interactions_validation >/dev/null 2>&1 ); ( cd $S/new && make -j2 crem_interactions crem_interactions_validation >/dev/null 2>&1 )
for v in old new; do [ -x $S/$v/crem_interactions ] && [ -x $S/$v/crem_interactions_validation ] || { echo "BLAD make $v"; exit 1; }; done
echo "build $(date +%T)"
for v in old new; do ( cd $S/$v && ./crem_interactions_validation > $D/out/validation_$v.txt 2>&1; echo EXIT $? >> $D/out/validation_$v.txt ); echo "validation $v $(date +%T)"; done
run(){ v=$1; tag=$2; shift 2; mkdir -p $S/reg/$v/$tag; t0=$(date +%s); ( cd $S/reg/$v/$tag && $S/${v%%_*}/crem_interactions "$@" > out.txt 2> err.txt ); echo "$v $tag $(( $(date +%s) - t0 ))" >> $D/out/walltime818.txt; }
: > $D/out/walltime818.txt
for v in old new new_off; do
  if [ $v = new_off ]; then export CREM_SEPARATION_FLOOR_SCALE=0.05; else unset CREM_SEPARATION_FLOOR_SCALE; fi
  run $v 1 --mode statistical --phenomenon 1 --runs 16 --seed 466
  run $v 2 --mode statistical --phenomenon 2 --runs 16 --seed 466
  run $v L2 --mode statistical --phenomenon 2 --runs 16 --seed 466 --level 2
  run $v 5Ps --mode statistical --phenomenon 5 --runs 50 --seed 551
  run $v 5H --mode statistical --phenomenon 5 --runs 50 --seed 551 --pair proton,electron
  run $v 6 --mode statistical --phenomenon 6 --runs 200 --seed 42
  echo "reg $v $(date +%T)"
done
unset CREM_SEPARATION_FLOOR_SCALE
for w in new new_off; do for tag in 1 2 L2 5Ps 5H 6; do s=SAME; for x in out.txt err.txt; do cmp -s $S/reg/old/$tag/$x $S/reg/$w/$tag/$x || s=DIFF; done; echo "old vs $w $tag $s"; done; done > $D/out/regression818.txt
mkdir -p $D/out/reg; for v in old new; do for tag in 1 2 L2 5Ps 6; do cp $S/reg/$v/$tag/out.txt $D/out/reg/${v}_$tag.out; done; done
echo "regression $(date +%T)"
# (4) Ps n = 2 orbit under the NEW default floor (no CREM_SEPARATION_FLOOR_SCALE), probe p814 from the new tree, setup of 815
O=$D/out/runs; mkdir -p $O/N $O/PN; export S O
B="MODE=est QSPIN=1 BUDGET=60 CREM_PAIR_TORQUE=200 CREM_PAIR_QM1S=1 CREM_DEBUG_BORN=1 CREM_EMISSION_LEDGER=1 CREM_PAIR_QM1S_TRACE=1e-12 BEXT=0.01 CREM_PAIR_TORQUE_ANYL=1 CREM_PAIR_TORQUE_APO=0.5 CREM_PAIR_TORQUE_PASS=1 CREM_SECULAR_NODES=1 NLEV=2 CREM_NO_DIPOLE_FORCE=1 CREM_PAIR_TORQUE_APSE=1 CREM_PAIR_TORQUE_PS=1 PAIR=electron,positron CREM_PAIR_TORQUE_CONSERVATIVE=1"
N5="u5a,0.5120091953982716 u5b,0.5652161499007791 u5c,0.66666666666666663 u5d,0.8125025652578265 u5e,0.9551918755393795"
{ for h in 4096 2048; do for nl in $N5; do IFS=, read nn f <<< "$nl"; echo "N/F1_ops_${nn}_s${h}_5 5 $B CREM_MECHANICAL_STEPS_PER_ORBIT=$h LFRAC=$f PHEN=2"; done; done
  for nl in $N5; do IFS=, read nn f <<< "$nl"; echo "PN/F0_pps_${nn}_s2048_5 5 $B CREM_MECHANICAL_STEPS_PER_ORBIT=2048 LFRAC=$f PHEN=1"; done
} | xargs -P 2 -L 1 bash -c 'ulimit -v 2500000; name=$0; seed=$1; shift
  env "$@" $S/p814_new $seed >/dev/null 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"'
echo "ALLDONE $(date +%T)"
