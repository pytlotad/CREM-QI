#!/bin/bash
cd $R652 && xargs -P 4 -I{} bash -c "{}" < jobs.txt
echo DONE
