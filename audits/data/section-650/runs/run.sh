#!/bin/bash
cd $R650 && xargs -P 4 -I{} bash -c "{}" < jobs.txt; echo DONE
