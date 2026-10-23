# Audit 578: capture distribution of the midpoint mode rule (577) for hydrogen at K/B = 0.05, n <= 8, and its exact QM cascade (557).
import math, numpy as np, mpmath as mp
exec(open('../section-557/cascade557.py').read().split('# capture shares at K/B = 0.05')[0])   # A, tot, ns, alpha, NMAX
src577=open("../section-577/pw577.py").read().rsplit("\nfor KoverB,label in",1)[0]
g={}; exec(src577,g)
KoverB=0.05; K=KoverB/2; kw=math.sqrt(2*K)
cap_rule={}; cap_qm={}
for n in range(1,NMAX+1):
    omega=kw*kw/2+1/(2*n*n)
    for lp in range(n):
        cap_qm[(n,lp)]=float(2*(2*lp+1)*(mp.mpf(omega)*ns['alpha']/mp.mpf(kw))**2*ns['sigma_PI'](n,lp,mp.mpf(kw)))
        s=0.0
        for lin in (lp-1,lp+1):
            if lin<0: continue
            s+=(2*lin+1)*g['rule_P'](lin,n,lp,K)
        cap_rule[(n,lp)]=s*math.pi/kw**2     # sigma = (pi/k^2) sum (2l+1) P
def cascade(cap):
    Z=sum(cap.values()); flux={k:v/Z for k,v in cap.items()}; lines={}
    for n in range(NMAX,1,-1):
        for l in range(n):
            f=flux.get((n,l),0)
            if f==0 or (n,l)==(2,0): continue
            for (a,b,m,lp),v in A.items():
                if (a,b)!=(n,l): continue
                p=f*v/tot[(n,l)]; flux[(m,lp)]=flux.get((m,lp),0)+p; lines[(n,m)]=lines.get((n,m),0)+p
    return dict(direct1s=cap[(1,0)]/Z, lya=lines.get((2,1),0), end2s=flux[(2,0)], ha=lines.get((3,2),0), hb=lines.get((4,2),0),
                nshare=[sum(v for (n,l),v in cap.items() if n==m)/Z for m in range(1,5)])
Zq=sum(cap_qm.values()); Zr=sum(cap_rule.values())
print(f"H K/B = {KoverB}: sigma(rule)/sigma(QM), n <= {NMAX}: {Zr/Zq:.3f}")
print("  (n, l): QM share / rule share / ratio")
for n in range(1,5):
    print("   "+'  '.join(f"{n}{'spdfghij'[l]} {cap_qm[(n,l)]/Zq:.4f}/{cap_rule[(n,l)]/Zr:.4f}/{(cap_rule[(n,l)]/Zr)/(cap_qm[(n,l)]/Zq):.2f}" for l in range(n)))
for name,cap in (('QM',cap_qm),('rule 577',cap_rule)):
    c=cascade(cap)
    print(f"  {name:9s}: direct 1s {c['direct1s']:.3f}, Lyman-alpha {c['lya']:.3f}, ends in 2s {c['end2s']:.3f}, Halpha {c['ha']:.4f}, Hbeta {c['hb']:.4f}, Halpha/Hbeta {c['ha']/c['hb']:.3f}; n shares "+' '.join(f"{x:.3f}" for x in c['nshare']))
print("  engine 553 (audit 557, split557): Halpha/Hbeta 1.75, ends in 2s 0.018; 554: 14, 0.127")
