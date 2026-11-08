#!/usr/bin/env python3
"""Audit 740, H740c: energy of the last 2p -> 1s photon of each run grouped by (j, F) of the 2p state and F' of the 1s state (ledger spin
directions, as ana739d.py); spread within each group and distances between groups.  usage: ana740c.py new.err"""
import sys, importlib.util, os
spec=importlib.util.spec_from_file_location('a',os.path.join(os.path.dirname(os.path.abspath(__file__)),'../section-739/ana739d.py'))
sys_argv=sys.argv; sys.argv=[sys.argv[0]]; a=importlib.util.module_from_spec(spec); spec.loader.exec_module(a); sys.argv=sys_argv
e=1.602176634e-19; h=6.62607015e-34; hc_nm=1239.8419843320026
def lines(path):
    out={}; idx=None
    for line in open(path):
        if line.startswith('START idx'): idx=int(line.split()[2])
        elif line.startswith('CREM_LEDGER t') and idx is not None:
            parts=line.split('|'); t=parts[0].split(); B=a.side(parts[1][3:]); A=a.side(parts[2][3:])
            if 1.0<=B['L']/a.hbar<2.0 and A['L']/a.hbar<1.0: out[idx]=(float(t[t.index('Eg')+1])/e,B,A)
    return out
g={}
for i,(E,B,A) in sorted(lines(sys.argv[1]).items()):
    j,F,c=a.basins(B,A); g.setdefault((j,F,0 if c<-0.5 else 1),[]).append(E)
print(f'{"group":24s} {"n":>3s} {"mean [eV]":>16s} {"spread [eV]":>12s} {"nm":>14s} {"MHz":>16s}')
for k in sorted(g):
    v=g[k]; m=sum(v)/len(v); s=max(v)-min(v)
    print(f'2P{"1/2" if k[0]==0.5 else "3/2"} F={k[1]} -> 1S F\'={k[2]}  {len(v):3d} {m:16.10f} {s:12.3e} {hc_nm/m:14.7f} {m*e/h/1e6:16.1f}')
