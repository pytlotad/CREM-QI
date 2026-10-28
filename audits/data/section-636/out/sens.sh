#!/bin/bash
cd /tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r636
run(){ tag=$1; shift; env "$@" ./p636_args; }
for v in def L635 L634; do
  case $v in def) E="";; L635) E="CREM_HIDDEN_J_ENDPOINT=1";; L634) E="CREM_HIDDEN_J_ENDPOINT=1 CREM_SPIN_LEGACY_L_NORM=1";; esac
  for i in 3 11 12; do echo "env $E NLEV=3 LFRAC=0.8333333333 ./p636 $i > sens_${v}_3d_$i.out 2>/dev/null"; done
  for i in 11 8 1; do echo "env $E NLEV=2 LFRAC=0.75 ./p636 $i > sens_${v}_2p_$i.out 2>/dev/null"; done
done | xargs -P 4 -I{} bash -c "{}"
echo DONE > sens_done.txt
