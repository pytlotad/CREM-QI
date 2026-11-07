#!/usr/bin/env python3
"""Audit 739, H739d: F' of the 1s state after the last 2p -> 1s photon of each run (rule Q from the ledger's spin directions), and (j, F) of
the 2p state before it.  usage: ana739d.py base.err new.err"""
import sys, math
import numpy as np
hbar=1.054571817e-34
QM={(0.5,0):0.0,(0.5,1):1/3,(1.5,1):2/3,(1.5,2):0.0}
def vec(t,key):
    i=t.index(key); return np.array(list(map(float,t[i+1:i+4])))
def side(s):
    t=s.split()
    return dict(E=float(t[t.index('E')+1]),L=float(t[t.index('L')+1]),ax=vec(t,'ax'),m1=vec(t,'m1'),m2=vec(t,'m2'))
def runs(path):
    out={}; idx=None
    for line in open(path):
        if line.startswith('START idx'): idx=int(line.split()[2])
        elif line.startswith('CREM_LEDGER t') and idx is not None:
            parts=line.split('|'); B=side(parts[1][3:]); A=side(parts[2][3:])
            if 1.0<=B['L']/hbar<2.0 and A['L']/hbar<1.0: out[idx]=(B,A)   # last p -> s photon
    return out
def basins(B,A):
    sp=B['m1']/np.linalg.norm(B['m1']); se=-B['m2']/np.linalg.norm(B['m2'])   # proton s || m, electron s || -m
    Lh=B['ax']/np.linalg.norm(B['ax']); l=1
    Jv=Lh*(l+0.5)+se; k=min(max(round(np.linalg.norm(Jv)),l),l+1); j=k-0.5
    Jh=Jv/np.linalg.norm(Jv); F=round(np.linalg.norm(Jh*k+sp)-0.5); F=min(max(F,round(abs(j-0.5))),round(j+0.5))
    ap=A['m1']/np.linalg.norm(A['m1']); ae=-A['m2']/np.linalg.norm(A['m2']); c=float(ap@ae)
    return j,F,c
for path in sys.argv[1:]:
    r=runs(path); n=len(r); n0=0; per={}
    for i,(B,A) in sorted(r.items()):
        j,F,c=basins(B,A); f0=c<-0.5; n0+=f0
        per.setdefault((j,F),[0,0,0]); per[(j,F)][0]+=1; per[(j,F)][1]+=f0; per[(j,F)][2]+=(c<0)
    p=n0/n if n else float('nan'); s=math.sqrt(0.25*0.75/n) if n else float('nan')
    print(f'{path}: {n} runs with a 2p -> 1s photon; F\'=0 (cos < -1/2) {n0}/{n} = {p:.4f}; QM 1/4, binomial sigma {s:.4f}, '
          f'deviation {(p-0.25)/s:+.2f} sigma')
    for key in sorted(per):
        m,k0,kc=per[key]
        print(f'   2P{"1/2" if key[0]==0.5 else "3/2"} F={key[1]}: {m:3d} runs, F\'=0 by rule Q {k0:3d} ({k0/m:.3f}; QM 6j {QM.get(key,float("nan")):.3f}),'
              f' cos < 0 (old singlet rule) {kc:3d} ({kc/m:.3f})')
