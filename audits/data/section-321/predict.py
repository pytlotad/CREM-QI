# First photon from the circular n=2 orbit (L = 2 hbar, E = -B/4), non-relativistic.
import numpy as np
rng=np.random.default_rng(7); N=4_000_000
c=rng.uniform(-1,1,N); keep=rng.uniform(0,2,N)<1+c*c; c=c[keep]
h=np.where(rng.uniform(0,1,len(c))<(1+c)**2/(2*(1+c*c)),1.0,-1.0)
Lp=np.sqrt(5-4*h*c)                       # |L - h hbar d| in hbar
ok=Lp<2; print("refused (L' >= 2, ceiling <= 0): %.4f"%(1-ok.mean()))
Lp=Lp[ok]; ceil=1/Lp**2-0.25              # in units of B
for name,Eg in (("hbar*omega",0.25),("Bohr 3/4",0.75)):
    E=np.minimum(Eg,ceil); trimmed=ceil<Eg
    Eafter=-0.25-E; nE=np.sqrt(-1/Eafter); e2=1-Lp**2/nE**2*1   # e^2 = 1 - L^2/n_E^2 (hbar units)
    print(f"{name:11s}: trimmed {trimmed.mean():.4f}; E_gamma/B mean {E.mean():.4f}; "
          f"E_gamma/(3B/4) mean {(E/0.75).mean():.4f}; n_E after: mean {nE.mean():.4f}, "
          f"quartiles {np.percentile(nE,25):.4f} {np.median(nE):.4f} {np.percentile(nE,75):.4f}; "
          f"P(n_E<1.05) {(nE<1.05).mean():.4f}; P(n_E=sqrt2 +-1e-9) {(abs(nE-2**0.5)<1e-9).mean():.4f}; <e^2> {e2.mean():.4f}")
print("L' after: quartiles %.4f %.4f %.4f, min %.4f"%(*np.percentile(Lp,[25,50,75]),Lp.min()))
