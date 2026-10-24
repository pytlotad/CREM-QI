#!/bin/bash
# H624c / H624d job: state chunk mode ; args: n l scale first last mode(new|legacy|old)
S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r624
cd /home/teddy/Projekty/CREM-QI
n=$1; l=$2; sc=$3; a=$4; b=$5; mode=$6
mkdir -p $S/C/$mode
bin=$S/p624; extra=""
if [ $mode = legacy ]; then export CREM_MECHANICAL_LEGACY_RATES=1 CREM_MECHANICAL_STATE_STREAM=1; fi
if [ $mode = old ]; then bin=$S/p624old; fi
PAIR=proton,electron NLEV=$n LHBAR=$l OBS=0.05 CREM_DEBUG_LADDER=1 CREM_MECHANICAL_HAZARD_SCALE=$sc \
  $bin $(seq $a $b) > $S/C/$mode/n${n}l${l}_${a}.out 2> $S/C/$mode/n${n}l${l}_${a}.err
