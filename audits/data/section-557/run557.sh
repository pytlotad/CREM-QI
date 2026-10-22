#!/bin/bash
# args: label index n L
label=$1; i=$2; n=$3; L=$4; S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r542/casc
mkdir -p $S/$label; cd $S/$label
if [ "$n" = "1" ]; then echo "SKIP n=1" > c$i.err; exit 0; fi
j=$(python3 -c "print(min(1.0,max(1e-3,$L/$n)))")
env CREM_EMISSION_LEDGER=1 PAIR=proton,electron NLEV=$n LFRAC=$j BUDGET=1800 $S/p505 $i 2> c$i.err
