#!/bin/bash
# Audit 844: static spin part (engine Z / ZD / ZDC vs BO'C c^-2) on the 838 node orbits scaled by s; o-Ps, p-Ps (positron spin flipped), H F0, H F1.
set -u; export LC_ALL=C; ulimit -v 4000000
R=/home/teddy/Projekty/CREM-QI; S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s844; D=$R/audits/data/section-844; O=$D/out/runs
rm -rf $S; mkdir -p $S $O
cd $R
g++ -O2 -std=c++20 $(root-config --cflags) -I. -Imodules $D/p844.cpp -o $S/p844 $(root-config --libs) 2>/dev/null || { echo "BLAD kompilacja"; exit 1; }
echo "build $(date +%T)"
while read sk nd pair rp e a b c x y z; do
  sys=$sk; [ $sk = ops ] && sys=ops
  for v in 1 2; do
    if [ $sk = ops ] && [ $v = 2 ]; then sys=pps; x=$(python3 -c "print(-float('$x'))"); y=$(python3 -c "print(-float('$y'))"); z=$(python3 -c "print(-float('$z'))")
    elif [ $v = 2 ]; then continue; fi
    for s in 1 2 4 8 16 32; do
      ( ulimit -v 2500000; $S/p844 $pair $rp $e $a $b $c $x $y $z 4096 $s ) > $O/${sys}_${nd}_s$s.txt || echo "FAIL $sys $nd $s"
    done
    for s in 1 32; do ( ulimit -v 2500000; $S/p844 $pair $rp $e $a $b $c $x $y $z 8192 $s ) > $O/${sys}_${nd}_s${s}_N2.txt || echo "FAIL N2 $sys $nd $s"; done
  done
done < $D/out/params839.txt
echo "ALLDONE $(date +%T)"
