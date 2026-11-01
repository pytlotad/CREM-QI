#!/bin/bash
# Audit 676: hypothesis checks on runs/ (run from this directory) -> out/check676.txt
D=..; e=1.602176634e-25
ext(){ grep -E 'PAIR eigen draw|PAIR field eigen \(1S\)' $1; grep -o 'CREM_LEDGER t [^ ]* Eg [^ ]*' $1; }
echo "H676a off/ = ../section-670/runs (bajt w bajt)"; n=0;m=0; for f in runs/off/*.err; do n=$((n+1)); cmp -s $f $D/section-670/runs/$(basename $f) && m=$((m+1)); done; echo "  $m/$n"
echo "H676b on/ |S1S*| per 1S state (ueV)"; cat runs/on/*.err | awk -v e=$e '/PAIR eigen draw/{for(i=1;i<=NF;i++)if($i=="P_S")p=$(i+1)} /PAIR eigen ledger/{for(i=1;i<=NF;i++){if($i=="S1S")s=$(i+1);if($i=="Zsa")z=$(i+1)} v=(s-z)/e; a=v<0?-v:v; if(a>mx)mx=a; k=sprintf("%.4f",p); n[k]++; x[k]=v; c++} END{for(k in n)printf "  P_S %s  %d photons  S1S* %.4f\n",k,n[k],x[k]; printf "  %d photons max %.4f  <= 3: %s\n",c,mx,(mx<=3?"TAK":"NIE")}'
echo "H676b on/ Eg + labels = ../section-662/runs"; n=0;m=0; for f in runs/on/*.err; do n=$((n+1)); [ "$(ext $f)" = "$(ext $D/section-662/runs/$(basename $f))" ] && m=$((m+1)); done; echo "  $m/$n"
echo "H676c b0/ import (ueV) and Eg + labels = ../section-670/b0"; cat runs/b0/*.err | awk -v e=$e '/PAIR eigen photon/{for(i=1;i<=NF;i++)if($i=="import")x=$(i+1)/e; a=x<0?-x:x; if(a>mx)mx=a; c++} END{printf "  %d photons max|import| %.5f\n",c,mx}'
n=0;m=0; for f in runs/b0/*.err; do n=$((n+1)); [ "$(ext $f)" = "$(ext $D/section-670/b0/$(basename $f))" ] && m=$((m+1)); done; echo "  $m/$n"
echo "H676d ps_on/ = ps_off/"; n=0;m=0; for f in runs/ps_on/*.err; do n=$((n+1)); cmp -s $f runs/ps_off/$(basename $f) && m=$((m+1)); done; echo "  $m/$n"
echo "opisowo: import on/ vs ../section-662/runs (ueV)"; for d in runs/on $D/section-662/runs; do cat $d/*.err | awk -v e=$e -v d=$d '/PAIR eigen photon/{for(i=1;i<=NF;i++)if($i=="import")x=$(i+1)/e; c++; s+=x*x; if(x<mn)mn=x; if(x>mx)mx=x} END{printf "  %s %d rms %.2f min %.2f max %.2f\n",d,c,sqrt(s/c),mn,mx}'; done
echo "opisowo: ziarna 0-7, on/ vs off/ (ueV)"; for d in on off; do cat $(ls runs/$d/A?_[0-7].err) | awk -v e=$e -v d=$d '/PAIR eigen ledger/{for(i=1;i<=NF;i++){if($i=="S1S")s=$(i+1);if($i=="Zsa")z=$(i+1);if($i=="Zsb")zb=$(i+1);if($i=="S2P")p=$(i+1)} x=(s-z)/e; y=(zb-z)/e; w=p/e; c++; q+=x*x; r+=y*y; t+=w*w} END{printf "  %s %d S1S* rms %.3f dZs rms %.3f S2P* rms %.3f\n",d,c,sqrt(q/c),sqrt(r/c),sqrt(t/c)}'; done
