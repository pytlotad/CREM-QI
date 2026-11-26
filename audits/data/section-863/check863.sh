#!/bin/bash
# Audit 863: byte-for-byte comparison of the 14 controls with the archives (830, 835, 836, 846, 862) and of production sets 2 / 5Ps with 861.
set -u; export LC_ALL=C; A=/home/teddy/Projekty/CREM-QI/audits/data; D=$A/section-863; X=$D/out/base; rm -rf $X; mkdir -p $X
S=/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/s863
for s in 830 835 836 846 862 861; do mkdir -p $X/$s; tar -xJf $A/section-$s/runs/runs$s.tar.xz -C $X/$s; done
n=0; r=0
cnt(){ grep -c "$1" $2; }
for pr in "N_O_ops 836/runs/OS/D_F1_ops_u5e_s2048" "N_O_pps 836/runs/OS/D_F0_pps_u5e_s2048" "N_O_H0 836/runs/OS/D_F0_H_u5e_s2048" "N_O_H1 836/runs/OS/D_F1_H_u5e_s2048" \
          "N_Z_ops 835/runs/ZS/D_F1_ops_u5e_s2048" "N_Z_pps 835/runs/ZS/D_F0_pps_u5e_s2048" "N_Z_H0 835/runs/ZS/D_F0_H_u5e_s2048" "N_Z_H1 835/runs/ZS/D_F1_H_u5e_s2048" \
          "V_O 846/runs/OZ/D_F1_ops_u5e_s2048" "V_Z 846/runs/ZDC/D_F1_ops_u5e_s2048" "L_O 830/runs/O/D_F1_ops_u5e_s2048" "L_Z 830/runs/Z/D_F1_ops_u5c_s2048" \
          "E_O 846/runs/OZ/D_F1_ops_u5e_s2048" "E_ZT 862/runs/ZDCT/D_F1_ops_u5e_s2048"; do set -- $pr
  if cmp -s $D/out/runs/$1.err $X/$2.err; then v=TAK; n=$((n+1)); else v=NIE; fi
  f=$D/out/runs/$1.err; l1=$(grep -m1 "PAIR START835" $f | awk '{print $4}')
  printf "%-8s = %-36s %s  (START835 %d lambda %s, SODIRAC841 %d, SOCANON843 %d, ORETZDC846 %d, CONS814 %d)\n" $1 $2 $v $(cnt START835 $f) "${l1:--}" $(cnt SODIRAC841 $f) $(cnt SOCANON843 $f) $(cnt ORETZDC846 $f) $(cnt CONS814 $f)
done
echo "H863k: $n/14 bajt w bajt -> $([ $n = 14 ] && echo TAK || echo NIE)"
for t in 2 5Ps; do same=TAK; for x in out.txt err.txt rc.txt; do cmp -s $S/reg/$t/$x $X/861/reg/new/$t/$x || same=NIE; done; echo "produkcja zestaw $t = 861 new: $same"; done
rm -rf $X
