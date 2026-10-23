# Audit 611: E1 J-branching -- photon as a rule-Q vector on the triangle cone (uniform azimuth) vs along L^ (610) vs QM 6j; scaling of the "against" branches.
import numpy as np, math
src=open('../section-610/j610.py').read().split('rng=np.random')[0]; g={}; exec(src,g); sixj=g['sixj']
print("T1  QM 'against' branch weights (spin 1):")
for fam,(dl,Jof,Jpof) in (("l->l-1, J=l-1 -> J'=l",(-1,lambda l:l-1,lambda l:l)),("l->l+1, J=l+1 -> J'=l",(+1,lambda l:l+1,lambda l:l))):
    vals=[]
    for l in (2,3,4,6,8,12,16,24,32,40):
        J=Jof(l); lp=l+dl
        w={jp:(2*jp+1)*sixj(l,J,1,jp,lp,1)**2 for jp in range(abs(lp-1),lp+2)}; z=sum(w.values()); p=w[Jpof(l)]/z
        vals.append((l,p))
    print(f"  {fam}: "+"  ".join(f"l={l}: {p:.2e} (l^2 p = {l*l*p:.3f}, l^4 p = {l**4*p:.2f})" for l,p in vals))
rng=np.random.default_rng(611)
def iso(n): v=rng.normal(size=(n,3)); return v/np.linalg.norm(v,axis=1,keepdims=True)
def branch(spin,l,dl,jph,N=400000):
    slen=spin+0.5; Lh=iso(N); sh=iso(N); L=(l+0.5)*Lh; s=slen*sh; lp=l+dl
    def lab(Lv,lv):
        x=np.linalg.norm(Lv+s,axis=1)-0.5; lo=abs(lv-spin); hi=lv+spin
        return np.clip(np.round(x-(spin%1))+(spin%1),lo,hi)
    J0=lab(L,l)
    if jph is None: Lp=(lp+0.5)*Lh
    else:
        Lm,Lpm=l+0.5,lp+0.5; c=(Lm*Lm+jph*jph-Lpm*Lpm)/(2*Lm*jph); c=np.clip(c,-1,1); sn=math.sqrt(1-c*c)
        t=np.where(np.abs(Lh[:,:1])<0.9,np.array([[1.,0,0]]),np.array([[0,1.,0]])); e1=np.cross(Lh,t); e1/=np.linalg.norm(e1,axis=1,keepdims=True); e2=np.cross(Lh,e1)
        ph=2*np.pi*rng.random(N)[:,None]; jv=jph*(c*Lh+sn*(np.cos(ph)*e1+np.sin(ph)*e2)); Lp=L-jv
    J1=lab(Lp,lp); out=[]
    for j in sorted(set(J0)):
        sel=J0==j; cl={jp:np.mean(J1[sel]==jp) for jp in np.arange(abs(lp-spin),lp+spin+0.1,1.0)}
        qm={jp:(2*jp+1)*sixj(l,j,spin,jp,lp,1)**2 for jp in cl}; z=sum(qm.values()); qm={k:v/z for k,v in qm.items()}
        out.append((j,cl,qm))
    return out
print("T2  max |P_klas - P_QM| per variant (along L^ = 610; triangle |j_ph| = 3/2; triangle sqrt 2):")
for tag,spin in (("H",0.5),("Ps",1.0)):
    for name,jph in (("along L^",None),("triangle 3/2",1.5),("triangle sqrt2",math.sqrt(2))):
        worst=0; worst_case=None
        for l in range(0 if tag=="H" else 1,5):
            for dl in (-1,1):
                if l+dl<0: continue
                for j,cl,qm in branch(spin,l,dl,jph):
                    for k in cl:
                        d=abs(cl[k]-qm[k])
                        if d>worst: worst=d; worst_case=(l,l+dl,j,k,cl[k],qm[k])
        print(f"  {tag:2s} {name:15s}: max {worst:.3f} at l={worst_case[0]}->{worst_case[1]} J={worst_case[2]} -> J'={worst_case[3]} ({worst_case[4]:.3f}/{worst_case[5]:.3f})",flush=True)
print("  examples (Ps, triangle 3/2): ",end="")
for j,cl,qm in branch(1.0,2,-1,1.5): print(f"3D{int(j)} -> "+", ".join(f"3P{int(k)} {cl[k]:.3f}/{qm[k]:.3f}" for k in cl),end="; ")
print()
