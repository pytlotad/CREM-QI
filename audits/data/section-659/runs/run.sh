#!/bin/bash
cd $R659 && xargs -P 4 -I{} bash -c "{}" < jobs.txt
echo DONE
