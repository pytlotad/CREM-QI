# Audit 549 (1b): only ALLOWED capture photons, K < hbar omega(r) <= K + B (n_cl >= 1) -- the stationary-phase reading of the rule (a photon hbar omega(r) leaves where the
# orbit's local frequency is omega); n = round(n_cl) (Kramers' band n +- 1/2).  Against Kramers sigma_n.  Hyperbola in the hyperbolic anomaly.
import numpy as np
eps0=8.8541878128e-12; e=1.602176634e-19; hbar=1.054571817e-34; c=299792458.0; me=9.1093837015e-31; mp=1.67262192369e-27; alpha=7.2973525693e-3
k=e*e/(4*np.pi*eps0); barn=1e-28
def run(name,mu,Ks,rcoll):
    a=hbar**2/(mu*k); B=k/(2*a)
    for KeV in Ks:
        K=KeV*e; v=np.sqrt(2*K/mu); ah=k/(2*K); T=np.sqrt(mu*ah**3/k)
        rlo=(hbar**2*k/(mu*(K+B)**2))**(1/3); rhi=(hbar**2*k/(mu*K*K))**(1/3)   # hbar omega = K+B at rlo, = K at rhi
        bs=np.geomspace(1e-14,2*rhi+np.sqrt(2*ah*2*rhi)*3,1500); db=np.gradient(bs); sig=0.0; W=[]; NC=[]
        for b,dbi in zip(bs,db):
            L=mu*v*b; ecc=np.sqrt(1+2*K*L*L/(mu*k*k)); rp=ah*(ecc-1)
            if rp<rcoll or rp>=rhi: continue
            Fm=np.arccosh((rhi/ah+1)/ecc); F=np.linspace(-Fm,Fm,6001); u=ecc*np.cosh(F)-1; r=ah*u
            om=np.sqrt(k/(mu*r**3)); P=e*e*(k/(mu*r*r))**2/(6*np.pi*eps0*c**3); dN=P/(hbar*om)*T*u
            ok=(r>=rlo)&(r<rhi); N=np.trapezoid(dN*ok,F); sig+=2*np.pi*b*N*dbi
            W.append(dN[ok]*(F[1]-F[0])*2*np.pi*b*dbi); NC.append(np.sqrt(B/(hbar*om[ok]-K)))
        w=np.concatenate(W); ncl=np.concatenate(NC); nr=np.maximum(1,np.rint(ncl)).astype(int)
        share=[w[nr==m].sum()/w.sum() for m in range(1,8)]
        eta2=B/K; ns=np.arange(1,20000); sk=(32*np.pi/(3*np.sqrt(3)))*alpha**3*a*a*eta2**2/(ns*(eta2+ns*ns))
        print(f"{name} K {KeV:5.2f} eV: allowed-photon sigma {sig/barn:.4e} b, Kramers {sk.sum()/barn:.4e} b, ratio {sig/sk.sum():.3f}; emission radii {rlo*1e12:.1f}-{rhi*1e12:.0f} pm")
        print(f"   n shares model: "+' '.join(f"{x:.3f}" for x in share)+" | Kramers: "+' '.join(f"{x:.3f}" for x in sk[:7]/sk.sum()))
        print(f"   per-n sigma ratio model/Kramers n=1..5: "+' '.join(f"{w[nr==m].sum()/sk[m-1]:.3f}" for m in range(1,6)))
run("Ps",me/2,(0.01,0.1,0.6,1.0),0.26e-12)
run("H ",me*mp/(me+mp),(0.01,0.1,0.6,1.0),0.26e-12)
