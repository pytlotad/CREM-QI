# Audit 602: is the classical 2P Ps spin-orbit-Zeeman system integrable? Largest Lyapunov exponent, 200 random initial conditions per field.
import numpy as np
A=2.7415; Bt=0.45725                     # GHz
mu=9.284765e-24; h=6.62607015e-34
TWO_PI=2*np.pi
def grads(L,s1,s2,a):
    Ln=np.linalg.norm(L,axis=1,keepdims=True); Lh=L/Ln
    S=s1+s2; d12=np.sum(s1*s2,1,keepdims=True); p1=np.sum(s1*Lh,1,keepdims=True); p2=np.sum(s2*Lh,1,keepdims=True)
    # H = A L.S + 4 Bt [ (3/2)(s1.s2 - p1 p2) - s1.s2 ] + 2a (s1z - s2z)   [GHz; time in ns -> rad via 2 pi]
    gs1=A*L+4*Bt*(1.5*(s2-p2*Lh)-s2); gs2=A*L+4*Bt*(1.5*(s1-p1*Lh)-s1)
    gs1[:,2]+=2*a; gs2[:,2]-=2*a
    # d(p1 p2)/dL = (p2 (s1 - p1 Lh) + p1 (s2 - p2 Lh))/|L|
    dpp=(p2*(s1-p1*Lh)+p1*(s2-p2*Lh))/Ln
    gL=A*S+4*Bt*(-1.5*dpp)
    return gL,gs1,gs2
def rhs(Y,a):
    L,s1,s2=Y[:,0:3],Y[:,3:6],Y[:,6:9]
    gL,g1,g2=grads(L,s1,s2,a)
    return TWO_PI*np.concatenate([np.cross(gL,L),np.cross(g1,s1),np.cross(g2,s2)],1)
def rk4(Y,a,dt):
    k1=rhs(Y,a); k2=rhs(Y+0.5*dt*k1,a); k3=rhs(Y+0.5*dt*k2,a); k4=rhs(Y+dt*k3,a)
    return Y+dt/6*(k1+2*k2+2*k3+k4)
def iso(rng,n):
    v=rng.normal(size=(n,3)); return v/np.linalg.norm(v,axis=1,keepdims=True)
def lyap(B,slen,n=200,T=100.0,dt=2e-3,renorm=0.5,seed=602):
    rng=np.random.default_rng(seed); a=mu*B/h/1e9
    Y=np.concatenate([1.5*iso(rng,n),slen*iso(rng,n),slen*iso(rng,n)],1)
    d0=1e-8; P=iso(rng,n*3).reshape(n,9); P/=np.linalg.norm(P,axis=1,keepdims=True); Z=Y+d0*P
    acc=np.zeros(n); steps=int(round(renorm/dt)); t=0.0
    while t<T-1e-9:
        for _ in range(steps): Y=rk4(Y,a,dt); Z=rk4(Z,a,dt)
        t+=renorm; D=Z-Y; d=np.linalg.norm(D,axis=1); acc+=np.log(d/d0); Z=Y+D*(d0/d)[:,None]
    return acc/T     # 1/ns
for slen,tag in ((1.0,'(i) |s| = 1 (rule Q)'),(0.5,'(ii) |s| = 1/2 (engine)')):
    base=lyap(0.0,slen); bmax=base.max()
    print(f"{tag}: baseline B = 0: lambda max {bmax:.4f} /ns, median {np.median(base):.4f}",flush=True)
    for B in (0.01,0.03,0.1,0.3,1.0,3.0):
        lam=lyap(B,slen); frac=np.mean(lam>5*bmax)
        print(f"   B {B:5.2f} T: lambda median {np.median(lam):.4f}, max {lam.max():.4f} /ns; chaotic (> 5 x baseline max) {100*frac:.1f} %",flush=True)
