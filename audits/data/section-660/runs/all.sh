#!/bin/bash
cd $R660/old && make -j2 crem_interactions >/dev/null 2>&1; cd $R660/new && make -j2 crem_interactions crem_interactions_validation >/dev/null 2>&1
cd $R660 && R=$R660 xargs -P 4 -I{} bash -c "{}" < jobs.txt
$R660/reg.sh
echo DONE
