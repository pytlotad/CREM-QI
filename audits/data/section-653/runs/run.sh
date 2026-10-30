#!/bin/bash
cd $R653 && xargs -P 4 -I{} bash -c "{}" < jobs.txt
echo DONE
