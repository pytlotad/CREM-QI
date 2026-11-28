#!/bin/bash
# Audit 874: prospective test of the w / r_peri rule (../section-871/local871.cpp unchanged); 15 runs, TEND = 6000 orbits.
set -u; export LC_ALL=C; ulimit -v 4000000
D=$(cd "$(dirname "$0")/.." && pwd); S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s874; B=$S/local871
mkdir -p $S; O=$D/out/runs; mkdir -p $O; g++ -O3 -march=native -std=c++17 -o $B $D/../section-871/local871.cpp || { echo BLAD kompilacja; exit 1; }
echo "build $(date +%T)"; export B O
PS="M1=1 G1=2.00231930436"
tend() { python3 -c "import math;m1=$1;mu=m1/(m1+1);print('%.0f'%(6000*2*math.pi*$2**3/mu))"; }
MH=1836.15267343
job() { local name=$1 m1=$2 n=$3 e=$4 w=$5; shift 5; local dx; dx=$(python3 -c "print('%.4f'%($w/5))")
  local sp=""; [ "$m1" = 1 ] && sp="$PS"; echo "$name $sp N0=$n E0=$e W=$w DX=$dx TEND=$(tend $m1 $n) $*"; }
{ job K        $MH 1.02 0.97 0.031 KAPPA=0
  job A_H97_q1 $MH 1.02 0.97 0.031
  job A_H97_q4 $MH 1.02 0.97 0.125
  job A_H85_q1 $MH 1.02 0.85 0.187
  job A_H85_q3 $MH 1.02 0.85 0.546
  job P_97_q1  1   1.02 0.97 0.062
  job P_97_q4  1   1.02 0.97 0.250
  job P_85_q1  1   1.02 0.85 0.375
  job P_85_q3  1   1.02 0.85 1.09
  job N2_H95   $MH 2.02 0.95 0.816
  job N2_H97   $MH 2.02 0.97 0.490
  job N2_P95   1   2.02 0.95 1.632
  job N3_H95   $MH 3.02 0.95 1.825
  job G1       $MH 1.02 0.95 0.25 GAM=8e-3
  job G2       $MH 2.02 0.95 0.816 GAM=1.25e-4; } > $S/jobs.txt
# longest first
python3 -c "
import sys;L=open('$S/jobs.txt').read().split('\n');L=[l.strip() for l in L if l.strip()]
L.sort(key=lambda l:-float(l.split('TEND=')[1].split()[0]));open('$S/jobs.txt','w').write('\n'.join(L)+'\n')"
echo "jobs $(wc -l < $S/jobs.txt)"; cat $S/jobs.txt
xargs -P 3 -L 1 bash -c 'ulimit -v 2500000; name=$0; env "$@" $B > $O/$name.txt 2>$O/$name.err && echo "done $name $(date +%T)" || echo "FAIL $name $(date +%T)"' < $S/jobs.txt
echo "ALLDONE $(date +%T)"
