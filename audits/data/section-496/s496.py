# Audit 496: rule Q (singlet iff cos th12 < -1/2, audit 476) vs Born (P_S = (1 - P+.P-)/4) for Ps formation on an optically pumped Rb
# target (Machacek et al., Atoms 11, 65 (2023) proposal).  Calibration-free observables: Sigma = P_S(par) + P_S(anti) - 1/2 (Born 0),
# Lambda = linearity of D = P_S(anti) - P_S(par) in P_Rb (Born 1).  Three classical decompositions of a partial polarization P.
# python3 s496.py (numpy).
import numpy as np
rng=np.random.default_rng(496); N=4_000_000
def iso(n):
    v=rng.normal(size=(n,3)); return v/np.linalg.norm(v,axis=1,keepdims=True)
def vmf(n,P):
    # von Mises-Fisher on S^2 about +z with mean cos = P: solve coth k - 1/k = P
    if P<=1e-9: return iso(n)
    if P>=1-1e-9: return np.tile([0,0,1.0],(n,1))
    lo,hi=1e-6,1e4
    for _ in range(200):
        k=0.5*(lo+hi); m=1/np.tanh(k)-1/k
        lo,hi=(k,hi) if m<P else (lo,k)
    u=rng.uniform(size=n); w=1+np.log(u+(1-u)*np.exp(-2*k))/k
    ph=rng.uniform(0,2*np.pi,n); s=np.sqrt(np.maximum(0,1-w*w))
    return np.stack([s*np.cos(ph),s*np.sin(ph),w],1)
def mix(n,P,cone=False):
    v=iso(n); m=rng.uniform(size=n)<P
    if cone:
        ph=rng.uniform(0,2*np.pi,m.sum()); c=0.5; s=np.sqrt(0.75)
        v[m]=np.stack([s*np.cos(ph),s*np.sin(ph),np.full(m.sum(),c)],1)
    else: v[m]=[0,0,1.0]
    return v
models={'(a) vMF':lambda P: vmf(N,P),'(b) axis+iso':lambda P: mix(N,P),'(c) cone60+iso':lambda P: mix(N,P,True)}
def ps_q(e,p,anti):
    if anti: p=p*np.array([1,1,-1.0])          # partner polarized along -z: mirror z
    return np.mean(np.sum(e*p,1)<-0.5)
print("check P+ = P_Rb = 0:",{k:round(ps_q(f(0),f(0),False),4) for k,f in models.items()})
PRb=[0,0.25,0.5,0.75,1.0]
for Pp in (0.3,0.6):
    print(f"\nP+ = {Pp}")
    print(f"  {'model':16s} " + " ".join(f"P_Rb={x:<4}: par / anti " for x in PRb))
    born=[((1-Pp*x)/4,(1+Pp*x)/4) for x in PRb]
    print(f"  {'Born':16s} " + " ".join(f"{a:.4f} / {b:.4f}      " for a,b in born))
    out={}
    for name,f in models.items():
        e=f(Pp); r=[]
        for x in PRb:
            p=f(x); r.append((ps_q(e,p,False),ps_q(e,p,True)))
        out[name]=r
        print(f"  {name:16s} " + " ".join(f"{a:.4f} / {b:.4f}      " for a,b in r))
    se=np.sqrt(0.25*0.75/N)
    print(f"  MC standard error per P_S ~ {se:.1e}")
    for name,r in [('Born',born)]+list(out.items()):
        D=[b-a for a,b in r]; S=[a+b-0.5 for a,b in r]
        Lam=(D[4]-D[0])/(2*(D[2]-D[0])) if abs(D[2]-D[0])>1e-12 else float('nan')
        print(f"  {name:16s} Sigma(P_Rb) = " + " ".join(f"{s:+.4f}" for s in S) + f"   Lambda {Lam:.3f}   D(1)/P+ {D[4]/Pp:.3f} (Born 0.5)")
    # statistical reach: project "3 % of P+ per hour" -> sigma(D) = 0.03 D_Born(1) per hour; sigma(Sigma) ~ sigma(D) (same counts)
    sD=0.03*Pp*0.5
    for name,r in out.items():
        S1=r[4][0]+r[4][1]-0.5
        h=(5*sD/abs(S1))**2 if abs(S1)>0 else float('inf')
        print(f"  {name:16s} |Sigma(1)| = {abs(S1):.4f} vs sigma ~ {sD:.4f}/sqrt(h): 5 sigma after {h:.2g} h")
