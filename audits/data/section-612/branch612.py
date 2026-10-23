# Audit 612 (report after the decomposition): branching ratios of np states -- engine (co-rotating harmonics, tail clipped to 1s; up from midpoint),
# 599 variant Go (geometric midpoint, both directions), and exact QM.
import numpy as np
exec(open('diag612.py').read().split('print("n   tau_eng')[0])
src599=open('../section-599/var599.py').read().split('rows=[]')[0]; g2={}; exec(src599,g2); rate599=g2['rate']
print("state  channel   QM share   engine share   Go share")
for n in (3,4,5,6):
    l=1; k,r=harmonics(n,l)
    ch=[(n2,0) for n2 in range(n-1,0,-1)]+[(n2,2) for n2 in range(n-1,2,-1)]
    qm={c:rate_qm(n,l,*c) for c in ch}
    eng={}
    for n2,l2 in ch:
        if l2==0: eng[(n2,l2)]=r[k==n-n2].sum() if n2>1 else r[k>=n-1].sum()
        else: eng[(n2,l2)]=rate_rule(n,l,n2,l2,'RBw')
    go={c:rate599('G','o',n,l,*c) for c in ch}
    tq,te,tg=sum(qm.values()),sum(eng.values()),sum(go.values())
    for c in ch:
        if qm[c]/tq<0.005: continue
        print(f"{n}p    -> {c[0]}{'sd'[c[1]//2]}     {qm[c]/tq:.3f}      {eng[c]/te:.3f}          {go[c]/tg:.3f}")
