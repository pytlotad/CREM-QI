# Audit 494 (author's T3): L on a cone -- does a wobbling orbit plane make 3P_1 = +2 reachable at QM lengths?
# sigma_1,2 = a S^ +- u, |u|^2 = sigma^2 - a^2, sigma^2 = 3; |L| = sqrt 2.  The latch u(p, psi) is optimized POINTWISE (2 deg scan) -- the
# most general u a single configuration can have.  Regime (1): orbit plane _|_ L(psi), Kepler azimuth uniform (r^-3 weight isotropic,
# 491a), wobble psi independent of the orbit.  Regime (2): r^ locked to the S^ azimuth (precession at the orbital frequency).
# python3 s494.py (numpy).
import numpy as np
rng=np.random.default_rng(494)
SIG2=3.0; LQ=np.sqrt(2.0)
NP=72; P=2*np.pi*np.arange(NP)/NP
AL=np.radians(np.arange(0,180,2.0))           # u and -u give the same S12, so 0..180 covers all
def frame(n):
    n=n/np.linalg.norm(n,axis=-1,keepdims=True)
    t=np.where(np.abs(n[...,2:3])<0.9,np.array([0,0,1.0]),np.array([1.0,0,0]))
    e1=np.cross(n,t); e1/=np.linalg.norm(e1,axis=-1,keepdims=True); return n,e1,np.cross(n,e1)
def s12_point(Sh,a,r):
    """max and min over the latch of 3(s1.r)(s2.r) - s1.s2; Sh, r (...,3); a (...)."""
    u=np.sqrt(np.maximum(SIG2-a*a,0.0))
    _,f1,f2=frame(Sh)
    U=(np.cos(AL)[:,None]*f1[...,None,:]+np.sin(AL)[:,None]*f2[...,None,:])*u[...,None,None]   # (...,nA,3)
    s1=a[...,None,None]*Sh[...,None,:]+U; s2=a[...,None,None]*Sh[...,None,:]-U
    v=3*np.sum(s1*r[...,None,:],-1)*np.sum(s2*r[...,None,:],-1)-np.sum(s1*s2,-1)
    return v.max(-1),v.min(-1)
def sph(th,ph): return np.stack([np.sin(th)*np.cos(ph),np.sin(th)*np.sin(ph),np.cos(th)],-1)
def config(NPSI=32,static=False):
    psi=2*np.pi*np.arange(NPSI)/NPSI
    k=lambda: rng.integers(1,4)
    thL=rng.uniform(0,np.pi); thS=rng.uniform(0,np.pi); chi=rng.uniform(0,2*np.pi)
    if static: dL=dS=fL=fS=g=0.0
    else: dL,dS,fL,fS=rng.uniform(0,1.2,4); g=rng.uniform(0,0.4) if rng.uniform()<0.5 else 0.0
    L=LQ*sph(thL+dL*np.sin(k()*psi+rng.uniform(0,6.3)),psi+fL*np.sin(k()*psi+rng.uniform(0,6.3)))
    Sh=sph(thS+dS*np.sin(k()*psi+rng.uniform(0,6.3)),psi+chi+fS*np.sin(k()*psi+rng.uniform(0,6.3)))
    a=np.sqrt(2.0)*(1+g*np.sin(psi+rng.uniform(0,6.3)))
    a=a*np.sqrt(2.0/np.mean(a*a)); a=np.minimum(a,np.sqrt(SIG2))            # <a^2> = 2 (S^2 = 2), a^2 <= sigma^2
    return L,Sh,a
def regime1(L,Sh,a):
    n,e1,e2=frame(L)
    r=np.cos(P)[None,:,None]*e1[:,None,:]+np.sin(P)[None,:,None]*e2[:,None,:]      # (psi, p, 3)
    hi,lo=s12_point(np.broadcast_to(Sh[:,None,:],r.shape),np.broadcast_to(a[:,None],r.shape[:2]),r)
    c=np.sum(n*Sh,-1); LS=np.mean(LQ*a*c)
    return hi.mean(),lo.mean(),LS,np.mean(c*c),np.mean(a*a*c*c),np.mean(a*a)
print("A. regime (1): wobble independent of the orbit; max/min <S12> over the pointwise latch")
rows=[]
for i in range(6000):
    L,Sh,a=config(static=(i<600)); hi,lo,LS,c2,ac2,a2=regime1(L,Sh,a)
    rows.append((LS,hi,lo,c2,ac2,i<600,a2))
R=np.array(rows)
bound=3-0.5*R[:,6]-1.5*R[:,4]          # pointwise-latch bound 3 - 0.5<a^2> - 1.5<a^2 c^2> with the sample's own <a^2>
viol=np.max(R[:,1]-bound)
full=R[:,6]>2-1e-9                     # samples that keep <a^2> = 2 (no clipping at sigma)
viol2=np.max(R[full,1]-(2-0.75*R[full,0]**2))
print(f"  {len(R)} configurations ({int(R[:,5].sum())} static); max (<S12>_max - [3 - 0.5<a^2> - 1.5<a^2 c^2>]) = {viol:+.4f};"
      f"  with <a^2> = 2 ({full.sum()}): max (<S12>_max - [2 - (3/4)<L.S>^2]) = {viol2:+.4f}")
for lsc in (-2.0,-1.0,0.0,1.0):
    m=np.abs(R[:,0]-lsc)<0.02
    if m.any():
        j=np.argmax(R[m,1]); print(f"  <L.S> = {lsc:+.0f} +- 0.02: {m.sum():4d} samples, envelope max <S12> {R[m,1].max():+.4f}"
              f" (sample <c^2> {R[m,3][j]:.3f}), min {R[m,2].min():+.4f};  analytic 2 - 3<L.S>^2/4 = {2-3*lsc*lsc/4:+.4f};"
              f" max over bin of (<S12>_max - own bound) {np.max(R[m,1]-bound[m]):+.4f}")
# exact constraint by construction: static configurations rotated to c = -1/2 with wobble added symmetric in c
print("  constrained scan at <L.S> = -1 exactly: c(psi) = -1/2 + w cos(psi), |S| fixed (wobble of the orbit plane about S)")
for w in (0.0,0.1,0.25,0.5):
    psi=2*np.pi*np.arange(64)/64; c=np.clip(-0.5+w*np.cos(psi),-1,1)
    Sh=np.tile([0,0,1.0],(64,1)); L=LQ*np.stack([np.sqrt(1-c*c)*np.cos(3*psi),np.sqrt(1-c*c)*np.sin(3*psi),c],-1)
    hi,lo,LS,c2,_,_=regime1(L,Sh,np.full(64,np.sqrt(2.0)))
    print(f"    w {w:.2f}: <L.S> {LS:+.4f}  <c^2> {c2:.4f}  max <S12> {hi:+.4f}  (2 - 3<c^2> = {2-3*c2:+.4f})  min {lo:+.4f}")
print("B. regime (2): r^ locked to the S^ azimuth (precession of L, S about J at the orbital frequency); c = -1/2")
c=-0.5; L=np.array([0,0,1.0]); Sh=np.array([np.sqrt(1-c*c),0,c])
best=[]
for p0 in np.radians(np.arange(0,181,5.0)):
    r=np.array([np.cos(p0),np.sin(p0),0.0])
    hi,lo=s12_point(Sh[None],np.array([np.sqrt(2.0)]),r[None]); best.append((np.degrees(p0),hi[0],lo[0]))
B=np.array(best)
print(f"  over the locked offset p0 (angle of r^ from the S^ projection in the orbit plane): max {B[:,1].max():+.3f} (p0 {B[np.argmax(B[:,1]),0]:.0f} deg),"
      f" min {B[:,2].min():+.3f};  max >= +2 for p0 in {B[B[:,1]>=2.0,0].min():.0f}..{B[(B[:,1]>=2.0)&(B[:,0]<=90),0].max():.0f} deg (and mirror)")
re=2.8179403262e-15; aPs=2*0.529177210903e-10; r2p=5*aPs   # <r> of 2p: (3n^2 - l(l+1))/2 a = 5 a_Ps
print(f"  needed precession rate / omega_orb ~ 1; the model's spin-orbit rate (3/4) r_e/r at <r>_2p = 5 a_Ps: {0.75*re/r2p:.2e}")
