#!/bin/bash
cd $R668/old && make -j2 crem_interactions >/dev/null 2>&1; cd $R668/new && make -j2 crem_interactions crem_interactions_validation >/dev/null 2>&1
cd $R668 && R=$R668 xargs -P 4 -I{} bash -c "{}" < jobs.txt
$R668/reg.sh
echo DONE
