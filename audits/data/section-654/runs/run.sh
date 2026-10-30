#!/bin/bash
cd $R654 && xargs -P 4 -I{} bash -c "{}" < jobs.txt
echo DONE
