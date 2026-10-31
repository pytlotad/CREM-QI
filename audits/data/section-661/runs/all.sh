#!/bin/bash
cd $R661/old && make -j2 crem_interactions >/dev/null 2>&1; cd $R661/new && make -j2 crem_interactions crem_interactions_validation >/dev/null 2>&1
cd $R661 && R=$R661 xargs -P 4 -I{} bash -c "{}" < jobs.txt
$R661/reg.sh
echo DONE
