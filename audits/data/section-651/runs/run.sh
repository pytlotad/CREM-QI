#!/bin/bash
cd $R651 && xargs -P 4 -I{} bash -c "{}" < jobs.txt; echo DONE
