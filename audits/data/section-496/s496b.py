# Audit 496 addendum: the observable in the apparatus of Machacek et al., Atoms 11, 65 (2023) -- delayed o-Ps counts only.
# S_o = (R+ + R-)/(2 R0) - 1 (R0: pump off; Born 0) and the Born-analysed polarimeter reading P_B,Born / P_B = 3 A / (P_T P_B).
# Same three classical decompositions as s496.py.  python3 s496b.py (numpy).
import numpy as np
rng=np.random.default_rng(4962); N=4_000_000
def iso(n):
    v=rng.normal(size=(n,3)); return v/np.linalg.norm(v,axis=1,keepdims=True)
def vmf(n,P):
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
        ph=rng.uniform(0,2*np.pi,m.sum()); s=np.sqrt(0.75)
        v[m]=np.stack([s*np.cos(ph),s*np.sin(ph),np.full(m.sum(),0.5)],1)
    else: v[m]=[0,0,1.0]
    return v
models={'(a) vMF':lambda P: vmf(N,P),'(b) axis+iso':lambda P: mix(N,P),'(c) cone60+iso':lambda P: mix(N,P,True)}
def ortho(e,p,anti):
    if anti: p=p*np.array([1,1,-1.0])
    return 1-np.mean(np.sum(e*p,1)<-0.5)
print("o-Ps fractions; R+ = target parallel to beam spin, R- = antiparallel, R0 = pump off")
for PB in (0.22,0.30):
    for PT in (0.9,1.0):
        print(f"\nP_B = {PB}, P_T = {PT}")
        rp=(3+PB*PT)/4; rm=(3-PB*PT)/4; r0=0.75
        print(f"  {'Born':15s} R+ {rp:.4f}  R- {rm:.4f}  R0 {r0:.4f}  S_o {(rp+rm)/(2*r0)-1:+.4f}  A {(rp-rm)/(rp+rm):.4f}  P_B(Born-analysed)/P_B {3*(rp-rm)/(rp+rm)/(PT*PB):.3f}")
        for name,f in models.items():
            e=f(PB); p=f(PT); z=f(0)
            rp=ortho(e,p,False); rm=ortho(e,p,True); r0=ortho(e,z,False)
            A=(rp-rm)/(rp+rm)
            print(f"  {name:15s} R+ {rp:.4f}  R- {rm:.4f}  R0 {r0:.4f}  S_o {(rp+rm)/(2*r0)-1:+.4f}  A {A:.4f}  P_B(Born-analysed)/P_B {3*A/(PT*PB):.3f}")
# statistics: 55 Hz delayed o-Ps signal (their estimate), T hours split equally over R+, R-, R0
for T in (1.0,):
    n=55*3600*T/3
    sS=np.sqrt(2/(4*n)+1/n)   # var of (R+ + R-)/(2 R0) - 1 for Poisson counts, equal expected counts
    print(f"\nstatistics: 55 Hz, {T:g} h split over three settings: sigma(S_o) ~ {sS:.4f} (Poisson only; beam normalization not included)")
