#!/bin/bash
# Audit 848: byte-for-byte comparison of the 7 controls with the archives (830, 835, 836, 846).
set -u; export LC_ALL=C; A=/home/teddy/Projekty/CREM-QI/audits/data; D=$A/section-848; X=$D/out/base; rm -rf $X; mkdir -p $X
for s in 830 835 836 846; do mkdir -p $X/$s; tar -xJf $A/section-$s/runs/runs$s.tar.xz -C $X/$s; done
n=0
for pr in "N_Z 846/runs/ZDC/D_F1_ops_u5e_s2048" "N_ZH 846/runs/ZDC/D_F0_H_u5e_s2048" "N_O 846/runs/OZ/D_F1_ops_u5e_s2048" \
          "L_Z 830/runs/Z/D_F1_ops_u5c_s2048" "L_O 830/runs/O/D_F1_ops_u5e_s2048" "T_ZS 835/runs/ZS/D_F1_ops_u5e_s2048" \
          "T_OS 836/runs/OS/D_F1_ops_u5e_s2048"; do set -- $pr
  if cmp -s $D/out/runs/$1.err $X/$2.err; then r=TAK; n=$((n+1)); else r=NIE; fi
  printf "%-5s = %-36s %s  (SODIRAC841 %d, SOCANON843 %d, ORETZDC846 %d)\n" $1 $2 $r $(grep -c SODIRAC841 $D/out/runs/$1.err) $(grep -c SOCANON843 $D/out/runs/$1.err) $(grep -c ORETZDC846 $D/out/runs/$1.err)
done
echo "H848k: $n/7 bajt w bajt -> $([ $n = 7 ] && echo TAK || echo NIE)"
rm -rf $X
