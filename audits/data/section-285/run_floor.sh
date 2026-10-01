#!/bin/bash
S=/tmp/claude-1000/-home-teddy-Projekty-Positronium/5981d9ee-dfbe-485a-8f93-5c6d87fd7653/scratchpad
IFS=: read -r floor steps <<< "$1"
TERSE=1 STEPS=$steps CREM_SEPARATION_FLOOR_SCALE=$floor $S/barrier2
