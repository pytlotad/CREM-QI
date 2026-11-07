#!/bin/bash
# usage: hazards.sh BIN OUTFILE -- first checkpoints' per-channel hazards (CREM_DEBUG_HAZARD_CHANNELS), wall budget 20 s per state
B=$1; F=$2; : > $F
for st in "H3d proton,electron 3 0.8333333333333334" "H4f proton,electron 4 0.875" "H4d proton,electron 4 0.625" "H4s proton,electron 4 0.125" \
          "H5s proton,electron 5 0.1" "H2p proton,electron 2 0.75" "Ps2p electron,positron 2 0.75" "Ps2s electron,positron 2 0.25"; do
  set -- $st; echo "== $1" >> $F
  env CREM_DEBUG_HAZARD_CHANNELS=1 PAIR=$2 NLEV=$3 LFRAC=$4 BUDGET=20 $B 0 2>&1 >/dev/null | grep HAZARD | sed -n '1,3p' >> $F
done
