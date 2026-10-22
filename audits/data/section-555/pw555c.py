# Audit 555 H2: other assignments of the partial wave l to the classical L (with the photon hbar), against exact QM.
import math, numpy as np
exec(open('pw555.py').read().split("for KoverB,label in")[0])
for KoverB,label in ((0.18,'Ps'),):
    K=KoverB/2; kw=math.sqrt(2*K); nmax=4
    q={lin:[qm_P(lin,n,kw) for n in range(1,nmax+1)] for lin in range(0,3)}
    for name in ('point l+1','band [l+1/2, l+3/2)'):
        print(f"{label} {name}:")
        for lin in range(0,3):
            if name.startswith('point'): c=cl_P(lin+1.0,K,nmax,True)[1:]
            else:
                Ls=np.linspace(lin+0.51,lin+1.49,25); c=sum(cl_P(L,K,nmax,True)*L for L in Ls)[1:]/Ls.sum()
            print(f"  l_in={lin}: model/QM n = 1..4: "+' '.join((f"{a/b:.2f}" if b>0 else ('inf' if a>0 else '-')) for a,b in zip(c,q[lin])))
