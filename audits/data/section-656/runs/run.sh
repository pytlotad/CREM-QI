#!/bin/bash
cd $R656 && xargs -P 4 -I{} bash -c "{}" < jobs.txt
echo DONE
