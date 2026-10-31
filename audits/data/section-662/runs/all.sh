#!/bin/bash
cd $R662/old && make -j2 crem_interactions >/dev/null 2>&1; cd $R662/new && make -j2 crem_interactions crem_interactions_validation >/dev/null 2>&1
cd $R662 && R=$R662 xargs -P 4 -I{} bash -c "{}" < jobs.txt
$R662/reg.sh
echo DONE
