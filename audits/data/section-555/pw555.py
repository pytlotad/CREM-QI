# Audit 555: partial-wave capture probabilities P(l_in -> n): exact QM (Gordon integrals) vs the model rule on the Coulomb hyperbola, with and without the photon hbar (audit 554).
# Atomic units of the pair (mu = hbar = k = 1): B = 1/2, K = k_w^2/2.
import sys, math, numpy as np, mpmath as mp
src=open('../section-553/gaunt553.py').read().split('print("check')[0]; ns={}; exec(src,ns)
alpha=float(ns['alpha'])
def qm_P(lin,n,kw):
    kw=mp.mpf(kw); omega=kw*kw/2+mp.mpf(1)/(2*n*n); tot=0
    for lp in (lin-1,lin+1):
        if lp<0 or lp>n-1: continue
        D=ns['dipole'](n,lp,kw,lin)
        sPI=4*mp.pi**2*ns['alpha']*omega/3*max(lp,lin)/mp.mpf(2*lp+1)*D**2
        srec=2*(2*lp+1)*(omega*ns['alpha']/kw)**2*sPI
        tot+=srec*kw*kw/(mp.pi*(2*lin+1))
    return float(tot)
def cl_P(L,K,nmax,spin):
    ah=1/(2*K); e=math.sqrt(1+2*K*L*L); F=np.linspace(-12,12,60001); u=e*np.cosh(F)-1; r=ah*u
    ok=r>0; F=F[ok]; u=u[ok]; r=r[ok]
    om=r**-1.5; P=(2/3)*alpha**3*r**-4; dN=P/om*np.sqrt(ah**3)*u; dF=F[1]-F[0]
    E=om; ncl=np.where(E>K,np.sqrt(0.5/np.maximum(E-K,1e-300)),0)
    out=np.zeros(nmax+1)
    for i in range(len(F)):
        if not (E[i]>K and ncl[i]>=0.5): continue
        n=max(1,int(round(ncl[i])))
        if n>nmax: continue
        if spin:
            s=-1.0 if L<1 else 1.0; Lt=abs(L-s)
            kin=0.5*(2*(K+1/r[i]))-(K+0.5/n/n)
            if not (Lt<n and kin>0 and Lt<=r[i]*math.sqrt(2*kin)): continue
        out[n]+=dN[i]*dF
    return out
for KoverB,label in ((0.18,'Ps'),(0.05,'H')):
    K=KoverB/2; kw=math.sqrt(2*K); nmax=4
    print(f"{label}: K/B = {KoverB} (k = {kw:.4f})")
    for lin in range(0,5):
        q=[qm_P(lin,n,kw) for n in range(1,nmax+1)]
        Ls=np.linspace(lin+0.01,lin+0.99,25); w=Ls
        c0=np.zeros(nmax+1); c1=np.zeros(nmax+1)
        for L in Ls: c0+=cl_P(L,K,nmax,False)*L; c1+=cl_P(L,K,nmax,True)*L
        c0/=w.sum(); c1/=w.sum()
        print(f"  l_in={lin}: QM P(n=1..4) "+' '.join(f"{x:.3e}" for x in q))
        print(f"          model no-hbar  "+' '.join(f"{x:.3e}" for x in c0[1:])+"   ratio "+' '.join((f"{a/b:.2f}" if b>0 else ('inf' if a>0 else '-')) for a,b in zip(c0[1:],q)))
        print(f"          model hbar     "+' '.join(f"{x:.3e}" for x in c1[1:])+"   ratio "+' '.join((f"{a/b:.2f}" if b>0 else ('inf' if a>0 else '-')) for a,b in zip(c1[1:],q)))
