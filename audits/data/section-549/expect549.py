# Audit 549: expected number of one-photon captures in experiment 5 for the sampling used (b ~ Rayleigh(sigma_b), K from the run's CAPDUMP lines), the rule of cap549c,
# hazard scale S: P(capture | b, K) = 1 - exp(-S N(b, K)).   python3 expect549.py CAPDUMPfile sigma_b_pm scale [H]
import sys, re, numpy as np
eps0=8.8541878128e-12; e=1.602176634e-19; hbar=1.054571817e-34; c=299792458.0; me=9.1093837015e-31
mp=1.67262192369e-27
k=e*e/(4*np.pi*eps0); mu=me*mp/(me+mp) if len(sys.argv)>4 and sys.argv[4]=='H' else me/2; a=hbar**2/(mu*k); B=k/(2*a)
def N_allowed(b,K):
    v=np.sqrt(2*K/mu); ah=k/(2*K); T=np.sqrt(mu*ah**3/k); rlo=(hbar**2*k/(mu*(K+4*B)**2))**(1/3); rhi=(hbar**2*k/(mu*K*K))**(1/3)
    L=mu*v*b; ecc=np.sqrt(1+2*K*L*L/(mu*k*k)); rp=ah*(ecc-1)
    if rp>=rhi: return 0.0
    Fm=np.arccosh((rhi/ah+1)/ecc); F=np.linspace(-Fm,Fm,4001); u=ecc*np.cosh(F)-1; r=ah*u
    om=np.sqrt(k/(mu*r**3)); P=e*e*(k/(mu*r*r))**2/(6*np.pi*eps0*c**3); dN=P/(hbar*om)*T*u
    return np.trapezoid(dN*(r>=rlo),F)
f,sb,S=sys.argv[1],float(sys.argv[2])*1e-12,float(sys.argv[3])
ev=[(float(re.search(r' K=(\S+)',l).group(1))*e,float(re.search(r' b=(\S+)',l).group(1))) for l in open(f) if l.startswith('CAPDUMP')]
p=[1-np.exp(-S*N_allowed(b,K)) for K,b in ev]
print(f"events {len(ev)}; expected one-photon captures {sum(p):.3e} (Poisson); max P per event {max(p):.3f}; mean S*N {np.mean([-np.log(1-x) for x in p]):.4f}")
