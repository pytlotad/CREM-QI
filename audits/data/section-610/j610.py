# Audit 610: E1 branching in J -- rule-Q basins (photon changes L by hbar along L^, spin fixed) vs QM 6j line strengths.
import numpy as np, math
def f(x): return math.factorial(int(round(x)))
def delta(a,b,c): return math.sqrt(f(a+b-c)*f(a-b+c)*f(-a+b+c)/f(a+b+c+1))
def tri(a,b,c): return abs(a-b)<=c<=a+b and abs(round(a+b+c)-(a+b+c))<1e-9
def sixj(a,b,c,d,e,F):
    if not (tri(a,b,c) and tri(a,e,F) and tri(d,b,F) and tri(d,e,c)): return 0.0
    pre=delta(a,b,c)*delta(a,e,F)*delta(d,b,F)*delta(d,e,c); s=0.0
    tmin=int(round(max(a+b+c,a+e+F,d+b+F,d+e+c))); tmax=int(round(min(a+b+d+e,b+c+e+F,c+a+F+d)))
    for t in range(tmin,tmax+1):
        s+=(-1)**t*f(t+1)/(f(t-a-b-c)*f(t-a-e-F)*f(t-d-b-F)*f(t-d-e-c)*f(a+b+d+e-t)*f(b+c+e+F-t)*f(c+a+F+d-t))
    return pre*s
rng=np.random.default_rng(610)
def iso(n): v=rng.normal(size=(n,3)); return v/np.linalg.norm(v,axis=1,keepdims=True)
def run(spin,l,dl,N=400000):
    slen=spin+0.5; Lh=iso(N); sh=iso(N); L=(l+0.5)*Lh; s=slen*sh
    def label(Lvec,lval):
        x=np.linalg.norm(Lvec+s,axis=1)-0.5          # J (or j) from |L + s| = (J + 1/2)
        lo=abs(lval-spin); hi=lval+spin
        return np.clip(np.round(x-(spin%1)) + (spin%1), lo, hi)    # half-integer j for spin 1/2, integer J for spin 1
    J0=label(L,l); lp=l+dl
    Lp=(lp+0.5)*Lh                                     # L' = L -+ hbar L^, along L^ (|L'| = l' + 1/2)
    J1=label(Lp,lp)
    out=[]
    for j in sorted(set(J0)):
        sel=J0==j; js=sorted(set(J1[sel]))
        cl={jp:np.mean(J1[sel]==jp) for jp in np.arange(abs(lp-spin),lp+spin+0.1,1.0)}
        qm={jp:(2*jp+1)*sixj(l,j,spin,jp,lp,1)**2 for jp in cl}; z=sum(qm.values()); qm={k:v/z for k,v in qm.items()}
        out.append((j,cl,qm))
    return out
worst={'H':0,'Ps':0}
for tag,spin in (('H',0.5),('Ps',1.0)):
    print(f"== {tag} (spin {spin}, length {spin+0.5})")
    for l in range(0 if tag=='H' else 1,5):
        for dl in (-1,+1):
            if l+dl<0 or (tag=='H' and l==0 and dl==-1): continue
            for j,cl,qm in run(spin,l,dl):
                d=max(abs(cl[k]-qm[k]) for k in cl); worst[tag]=max(worst[tag],d)
                print(f"  l={l} -> {l+dl}, j={j}: "+"  ".join(f"j'={k}: {cl[k]:.3f}/{qm[k]:.3f}" for k in cl)+f"   max diff {d:.3f}")
print(f"max diff H {worst['H']:.3f}, Ps {worst['Ps']:.3f}")
