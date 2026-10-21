# Audit 548: captured level n_eff = sqrt(Ry/|E_Coulomb,end|) and L/hbar in experiment 5 (CAPDUMP), against Kramers' level shares at each captured event's K
# (sigma_n ~ eta^4/(n (eta^2 + n^2)), eta^2 = Ry/K; Stobbe exceeds Kramers by 20/12/7.5 % at n = 1/2/4 -- Kotelnikov, Milstein 2019).  python3 ana548.py file Ry_eV label
import math, re, sys
f, Ry, label = sys.argv[1], float(sys.argv[2]), sys.argv[3]
cap=[]; outs={}
for line in open(f):
    if not line.startswith('CAPDUMP'): continue
    out=re.search(r'out=(\S+)',line).group(1); outs[out]=outs.get(out,0)+1
    if out in ('Scattering','Collision') or 'nresolved' in out: continue
    K=float(re.search(r' K=(\S+)',line).group(1)); E=float(re.search(r'Ecoul=\S+/(\S+)',line).group(1)); L=float(re.search(r' L=(\S+)',line).group(1))
    if E<0: cap.append((K,math.sqrt(Ry/abs(E)),L))
print(f"{label}: outcomes {outs}; captured (E_end < 0) {len(cap)}")
if not cap: sys.exit()
nmax=12; hist=[0]*(nmax+1); kram=[0.0]*(nmax+1)
for K,n,L in cap:
    hist[min(nmax,max(1,int(round(n))))]+=1
    eta2=Ry/K; w=[eta2**2/(m*(eta2+m*m)) for m in range(1,2000)]; s=sum(w)
    for m in range(1,nmax): kram[m]+=w[m-1]/s
    kram[nmax]+=sum(w[nmax-1:])/s
N=len(cap)
print("  n (rounded n_eff):    "+' '.join(f"{m:>5d}" for m in range(1,nmax))+f"  >={nmax}")
print("  model share:          "+' '.join(f"{hist[m]/N:5.3f}" for m in range(1,nmax+1)))
print("  Kramers share:        "+' '.join(f"{kram[m]/N:5.3f}" for m in range(1,nmax+1)))
ns=sorted(c[1] for c in cap); Ls=sorted(c[2] for c in cap); Ks=sorted(c[0] for c in cap)
med=lambda v: v[len(v)//2]
print(f"  n_eff min/median/max {ns[0]:.2f}/{med(ns):.2f}/{ns[-1]:.2f}; L/hbar min/median/max {Ls[0]:.3f}/{med(Ls):.3f}/{Ls[-1]:.3f}; K_CM of captures median {med(Ks):.3f} eV")
