# Audit 549 (1): the stochastic-mode photon rule on an unperturbed Coulomb hyperbola -- photons at rate P/(hbar omega), energy hbar omega(r), omega = sqrt(k/(mu r^3));
# binding photons (hbar omega > K): cross section and level n_cl = sqrt(B/(hbar omega - K)), rounded; against Kramers sigma_n = (32 pi/3 sqrt3) alpha^3 a^2 eta^4/(n (eta^2 + n^2)).
# Hyperbola in the hyperbolic anomaly: r = a_h (e cosh F - 1), dt = sqrt(mu a_h^3/k) (e cosh F - 1) dF, a_h = k/(2K).
import numpy as np
eps0=8.8541878128e-12; e=1.602176634e-19; hbar=1.054571817e-34; c=299792458.0; me=9.1093837015e-31; mp=1.67262192369e-27; alpha=7.2973525693e-3
k=e*e/(4*np.pi*eps0); barn=1e-28
def run(name,mu,Ks):
    a=hbar**2/(mu*k); B=k/(2*a)
    for KeV in Ks:
        K=KeV*e; v=np.sqrt(2*K/mu); rK=(hbar**2*k/(mu*K*K))**(1/3); ah=k/(2*K); T=np.sqrt(mu*ah**3/k)
        bs=np.geomspace(1e-16,1e-10,800); db=np.gradient(bs); sig=0.0; allw=[]; alln=[]; Nmax=0
        for b,dbi in zip(bs,db):
            L=mu*v*b; ecc=np.sqrt(1+2*K*L*L/(mu*k*k))
            x=rK/ah                                  # e cosh F - 1 < x
            if ecc-1>=x: continue
            Fm=np.arccosh((x+1)/ecc); F=np.linspace(-Fm,Fm,3001); u=ecc*np.cosh(F)-1; r=ah*u
            om=np.sqrt(k/(mu*r**3)); P=e*e*(k/(mu*r*r))**2/(6*np.pi*eps0*c**3)
            dN=P/(hbar*om)*T*u; N=np.trapezoid(dN,F); Nmax=max(Nmax,N)
            sig+=2*np.pi*b*N*dbi
            w=dN*(F[1]-F[0])*2*np.pi*b*dbi; ncl=np.sqrt(B/np.maximum(hbar*om-K,1e-300))
            allw.append(w); alln.append(ncl)
        w=np.concatenate(allw); ncl=np.concatenate(alln)
        nr=np.maximum(1,np.rint(ncl)).astype(int); share=[w[nr==m].sum()/w.sum() for m in range(1,8)]
        eta2=B/K; ns=np.arange(1,20000); sk=(32*np.pi/(3*np.sqrt(3)))*alpha**3*a*a*eta2**2/(ns*(eta2+ns*ns)); kshare=sk[:7]/sk.sum()
        print(f"{name} K {KeV:5.2f} eV: model-rule sigma {sig/barn:.4e} b, Kramers {sk.sum()/barn:.4e} b, ratio {sig/sk.sum():.3f}; max photons per pass {Nmax:.2e}; n_cl < 1 share {w[ncl<0.5].sum()/w.sum():.3f}")
        print(f"   n shares (rounded n_cl, n<1 -> 1) model: "+' '.join(f"{x:.3f}" for x in share)+" | Kramers: "+' '.join(f"{x:.3f}" for x in kshare))
run("Ps",me/2,(0.01,0.1,0.6,1.0))
run("H ",me*mp/(me+mp),(0.01,0.1,0.6,1.0))
