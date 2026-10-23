# Audit 579: expected (Rao-Blackwellized) channel distribution per forced event -- the rule's P(n, l') recomputed offline (pw577.rule_P,
# = C++ to 3e-4) from K and b of each FORCEDEVENT (L = mu v b), weighted HT x P_channel; against the drawn channels (FORCED5); QM cascade (H).
import sys, math, os, io, contextlib, numpy as np
src=open('../section-577/pw577.py').read().rsplit("\nfor KoverB,label in",1)[0]
os.chdir('../section-577'); g={}; exec(src,g); os.chdir('../section-579')
e=1.602176634e-19; hbar=1.054571817e-34; eps0=8.8541878128e-12; me=9.1093837015e-31; mp=1.67262192369e-27; kC=e*e/(4*math.pi*eps0)
errf,pair=sys.argv[1],sys.argv[2]
mu=me/2 if pair=='Ps' else me*mp/(me+mp); Eunit=mu*kC*kC/hbar**2
def sigb(K):
    Kj=K*e; rK=(hbar**2*kC/(mu*Kj*Kj))**(1/3); return 0.5*math.sqrt(rK*(rK+kC/Kj))
exp={}; drawn={}; nev=0
lines=open(errf).read().splitlines()
f5=[l for l in lines if l.startswith('FORCED5') and 'mode-rule' in l]
for l in lines:
    if not l.startswith('FORCEDEVENT'): continue
    d=dict(x.split('=',1) for x in l.split()[1:]); K=float(d['K']); b=float(d['b'])
    v=math.sqrt(2*K*e/mu); L=mu*v*b/hbar; lin=int(math.floor(L)); Kau=K*e/Eunit
    s=sigb(K); HT=2*math.pi*s*s*math.exp(b*b/(2*s*s))
    for lp in (lin-1,lin+1):
        if lp<0 or lp+1>16: continue
        for n in range(lp+1,17):
            P=g['rule_P'](lin,n,lp,Kau)
            exp[(n,lp)]=exp.get((n,lp),0)+HT*P
    nev+=1
for l in f5:
    d=dict(x.split('=',1) for x in l.split()[1:] if '=' in x); key=(int(d['n']),int(d['l'])); drawn[key]=drawn.get(key,0)+1
Z=sum(exp.values()); Zd=sum(drawn.values())
print(f"{pair}: {nev} events; expected (HT x P_channel) vs drawn-count shares (unweighted draws, conditional on each event):")
for n in range(1,5):
    print("  "+'  '.join(f"{n}{'spdfg'[l]} {exp.get((n,l),0)/Z:.4f}" for l in range(n)))
print(f"  n > 4 expected {sum(v for (n,l),v in exp.items() if n>4)/Z:.4f}")
if pair=='H':
    os.chdir('../section-557'); sc=open('cascade557.py').read().split('# capture shares at K/B = 0.05')[0]; h={}
    with contextlib.redirect_stdout(io.StringIO()): exec(sc,h)
    os.chdir('../section-579'); A,tot,NMAX=h['A'],h['tot'],h['NMAX']
    cap={kk:v for kk,v in exp.items() if kk[0]<=NMAX}; Zc=sum(cap.values()); flux={kk:v/Zc for kk,v in cap.items()}; ln={}
    for n in range(NMAX,1,-1):
        for l in range(n):
            f=flux.get((n,l),0)
            if f==0 or (n,l)==(2,0): continue
            for (a,b_,m,lp),v in A.items():
                if (a,b_)!=(n,l): continue
                p=f*v/tot[(n,l)]; flux[(m,lp)]=flux.get((m,lp),0)+p; ln[(n,m)]=ln.get((n,m),0)+p
    print(f"  QM cascade from the expected distribution of these events (n <= {NMAX}): direct 1s {cap.get((1,0),0)/Zc:.3f}, Lyman-alpha {ln.get((2,1),0):.3f}, ends in 2s {flux.get((2,0),0):.3f}, Halpha/Hbeta {ln.get((3,2),0)/ln.get((4,2),1e-300):.3f}")
