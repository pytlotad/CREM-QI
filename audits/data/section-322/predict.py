# First photon from n = 2 energy with L = hbar (e^2 = 3/4), non-relativistic, units B, hbar.
import numpy as np
rng=np.random.default_rng(11); N=4_000_000
c=rng.uniform(-1,1,N); c=c[rng.uniform(0,2,N)<1+c*c]
h=np.where(rng.uniform(0,1,len(c))<(1+c)**2/(2*(1+c*c)),1.0,-1.0)
Lp2=2*(1-h*c)                                   # |L - h d|^2 with |L| = 1
ok=Lp2<4; print("refused (L' >= 2): %.5f"%(1-ok.mean()))
Lp2=Lp2[ok]; ceil=1/Lp2-0.25
print("L' quartiles %.4f %.4f %.4f; P(L' <= 1) = P(hc >= 1/2) = %.4f"%(*np.sqrt(np.percentile(Lp2,[25,50,75])),(Lp2<=1).mean()))
for name,Eg in [("Bohr 3/4",0.75)]+[(f"hbar*omega, k={k}",0.25*k) for k in (1,2,3,4,5)]:
    E=np.minimum(Eg,ceil); tr=ceil<Eg; nE=np.sqrt(1/(0.25+E))
    print(f"{name:16s}: trimmed {tr.mean():.4f}; n_E after median {np.median(nE):.4f}; "
          f"P(n_E <= 1+1e-6) {(nE<=1+1e-6).mean():.4f}; P(n_E < 1) {(nE<1-1e-6).mean():.4f}; E/(3B/4) mean {(E/0.75).mean():.4f}")
