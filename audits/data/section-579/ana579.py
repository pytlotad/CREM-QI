# Audit 579: (n, l) distribution of the forced captures (FORCEDEVENT, weights HT x P as in ana551), the s -> 1s count (FORCED5 lines),
# and the exact QM cascade (A from cascade557) from the engine's distribution (hydrogen).  python3 ana579.py out.txt err.txt B_eV pair
import sys, re, math, numpy as np
out,err,B,pair=sys.argv[1],sys.argv[2],float(sys.argv[3]),sys.argv[4]
e=1.602176634e-19; hbar=1.054571817e-34; eps0=8.8541878128e-12; me=9.1093837015e-31; mp=1.67262192369e-27
k=e*e/(4*math.pi*eps0); mu=me/2 if pair=='Ps' else me*mp/(me+mp)
def sigb(K):
    Kj=K*e; rK=(hbar**2*k/(mu*Kj*Kj))**(1/3); return 0.5*math.sqrt(rK*(rK+k/Kj))
HT=lambda b,K: 2*math.pi*sigb(K)**2*math.exp(b*b/(2*sigb(K)**2))
W={}; s1s=0; nF=0
for l in open(err):
    if l.startswith('FORCED5') and 'mode-rule' in l:
        nF+=1; d=dict(x.split('=',1) for x in l.split()[1:] if '=' in x)
        if int(d['l_in'])==0 and int(d['n'])==1: s1s+=1
    if l.startswith('FORCEDEVENT'):
        d=dict(x.split('=',1) for x in l.split()[1:])
        if d['captured'] in ('Para-Positronium','Ortho-Positronium') and float(d['E'])<0:
            n=max(1,round(math.sqrt(B/abs(float(d['E']))))); L=float(d['L']); l_=max(0,min(n-1,round(L-0.5)))
            w=float(d['P'])*HT(float(d['b']),float(d['K'])); W[(n,l_)]=W.get((n,l_),0)+w
Z=sum(W.values())
print(f"{pair}: forced photons (mode rule) {nF}, s -> 1s {s1s}; weighted (n, l) shares:")
for n in range(1,6):
    print("  "+'  '.join(f"{n}{'spdfgh'[l]} {W.get((n,l),0)/Z:.4f}" for l in range(n)))
print(f"  n > 5: {sum(v for (n,l),v in W.items() if n>5)/Z:.4f}")
if pair=='H':
    import os; here=os.getcwd(); os.chdir('../section-557')
    src=open('cascade557.py').read().split('# capture shares at K/B = 0.05')[0]; g={}; 
    import io, contextlib
    with contextlib.redirect_stdout(io.StringIO()): exec(src,g)
    os.chdir(here)
    A,tot,NMAX=g['A'],g['tot'],g['NMAX']
    cap={kk:v for kk,v in W.items() if kk[0]<=NMAX}; Zc=sum(cap.values()); flux={kk:v/Zc for kk,v in cap.items()}; lines={}
    for n in range(NMAX,1,-1):
        for l in range(n):
            f=flux.get((n,l),0)
            if f==0 or (n,l)==(2,0): continue
            for (a,b,m,lp),v in A.items():
                if (a,b)!=(n,l): continue
                p=f*v/tot[(n,l)]; flux[(m,lp)]=flux.get((m,lp),0)+p; lines[(n,m)]=lines.get((n,m),0)+p
    ha,hb=lines.get((3,2),0),lines.get((4,2),0)
    print(f"  QM cascade from the engine's distribution (n <= {NMAX}): direct 1s {cap.get((1,0),0)/Zc:.3f}, Lyman-alpha {lines.get((2,1),0):.3f}, ends in 2s {flux.get((2,0),0):.3f}, Halpha/Hbeta {ha/hb if hb>0 else float('nan'):.3f}  (QM at K/B 0.05: 0.447 / 0.332 / 0.085 / 2.985; engine 553: 2s 0.018, Ha/Hb 1.75)")
