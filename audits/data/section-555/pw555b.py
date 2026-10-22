# Audit 555 H1: classical capture probability evaluated at the Langer point L = (l + 1/2) hbar of each partial wave, with the photon hbar (554), against exact QM.
import math, numpy as np
exec(open('pw555.py').read().split("for KoverB,label in")[0])
for KoverB,label in ((0.18,'Ps'),(0.05,'H')):
    K=KoverB/2; kw=math.sqrt(2*K); nmax=4
    print(f"{label}: K/B = {KoverB}")
    for lin in range(0,4):
        q=[qm_P(lin,n,kw) for n in range(1,nmax+1)]
        c=cl_P(lin+0.5,K,nmax,True)[1:]
        print(f"  l_in={lin}: Langer-point model / QM, n = 1..4: "+' '.join((f"{a/b:.2f}" if b>0 else ('inf' if a>0 else '-')) for a,b in zip(c,q)))
