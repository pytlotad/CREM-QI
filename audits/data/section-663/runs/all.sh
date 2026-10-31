#!/bin/bash
cd $R663/old && make -j2 crem_interactions >/dev/null 2>&1; cd $R663/new && make -j2 crem_interactions crem_interactions_validation >/dev/null 2>&1
cd $R663 && R=$R663 xargs -P 4 -I{} bash -c "{}" < jobs.txt
$R663/reg.sh
echo DONE
