#!/bin/bash
# audit 742 census (print only, CREM_LEVEL_WINDOW_CENSUS=1): p585 probes, p623 hydrogen i = 0..7, statistical scenarios of reg.sh
S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s742; T=$S/tree
rm -rf $T; mkdir -p $T; cd /home/teddy/Projekty/CREM-QI && rsync -a --exclude .git --exclude audits ./ $T/ && cp -r audits/data/section-585 audits/data/section-623 $T/
( cd /home/teddy/Projekty/CREM-QI && for p in 585 623; do g++ -O2 -std=c++20 $(root-config --cflags) -I. -Imodules audits/data/section-$p/p$p.cpp -o $S/p$p $(root-config --libs); done )
( cd $T && make -j4 crem_interactions >/dev/null 2>&1 )
export CREM_LEVEL_WINDOW_CENSUS=1
mkdir -p $S/probes; cd $S/probes
for spec in "H2p PAIR=proton,electron NLEV=2 LFRAC=0.75" "H3p PAIR=proton,electron NLEV=3 LFRAC=0.5" "H4p PAIR=proton,electron NLEV=4 LFRAC=0.375" \
  "H3s PAIR=proton,electron NLEV=3 LFRAC=0.16666666666666666" "H3d PAIR=proton,electron NLEV=3 LFRAC=0.8333333333333334" "Ps2p NLEV=2 LFRAC=0.75"; do
  set -- $spec; tag=$1; shift; env CREM_EMISSION_LEDGER=1 "$@" $S/p585 0 > $tag.out 2> $tag.err; echo "probe $tag $(grep LEVEL_WINDOW742 $tag.err)"
done > $S/census.txt
( mkdir -p $S/p623 && cd $S/p623 && env MODE=est PAIR=proton,electron PHEN=1 QSPIN=1 NLEV=2 LFRAC=0.75 BUDGET=3000 $S/p623 $(seq 0 7) > o.txt 2> e.txt; echo "p623 i=0..7 $(grep LEVEL_WINDOW742 e.txt)" ) >> $S/census.txt
run(){ tag=$1; shift; mkdir -p $S/reg/$tag; ( cd $S/reg/$tag && $T/crem_interactions "$@" > out.txt 2> err.txt ); echo "reg $tag $(grep LEVEL_WINDOW742 $S/reg/$tag/err.txt)" >> $S/census.txt; }
run 1 --mode statistical --phenomenon 1 --runs 16 --seed 466
run 2 --mode statistical --phenomenon 2 --runs 16 --seed 466
run L2 --mode statistical --phenomenon 2 --runs 16 --seed 466 --level 2
run 5Ps --mode statistical --phenomenon 5 --runs 50 --seed 551
run 5H --mode statistical --phenomenon 5 --runs 50 --seed 551 --pair proton,electron
run 6 --mode statistical --phenomenon 6 --runs 200 --seed 42
echo ALLDONE >> $S/census.txt
