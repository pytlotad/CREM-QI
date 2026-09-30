#!/bin/bash
S=/tmp/claude-1000/-home-teddy-Projekty-Positronium/5981d9ee-dfbe-485a-8f93-5c6d87fd7653/scratchpad
sd=$1
echo "OFF $($S/e2e_p5 $sd 120 1 | head -1)"
echo "ON  $(CREM_COMPUTED_EMISSION_PATTERN=1 $S/e2e_p5 $sd 120 1 | head -1)"
