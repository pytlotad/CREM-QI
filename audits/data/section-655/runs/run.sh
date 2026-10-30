#!/bin/bash
cd $R655 && xargs -P 4 -I{} bash -c "{}" < jobs.txt
echo DONE
