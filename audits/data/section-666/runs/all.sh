#!/bin/bash
cd $R666/old && make -j2 crem_interactions >/dev/null 2>&1; cd $R666/new && make -j2 crem_interactions crem_interactions_validation >/dev/null 2>&1
cd $R666 && R=$R666 xargs -P 4 -I{} bash -c "{}" < jobs.txt
$R666/reg.sh
echo DONE
