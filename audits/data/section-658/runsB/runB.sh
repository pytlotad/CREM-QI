#!/bin/bash
cd $R658 && xargs -P 4 -I{} bash -c "{}" < jobsB.txt
echo DONE
