#!/bin/bash
cd $R657 && xargs -P 4 -I{} bash -c "{}" < jobs.txt
echo DONE
