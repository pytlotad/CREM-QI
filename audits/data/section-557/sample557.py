# Audit 557: resample captured states (n, L) from the forced-capture lines with weights HT x P / sigma_Kramers (as ana551/lmodel554), n <= 8; writes "index n L" lines.
import sys, re, math, random, numpy as np
err,B,label,M=sys.argv[1],float(sys.argv[2]),sys.argv[3],int(sys.argv[4])
e=1.602176634e-19; hbar=1.054571817e-34; me=9.1093837015e-31; mp_=1.67262192369e-27; eps0=8.8541878128e-12; alpha=7.2973525693e-3
k=e*e/(4*np.pi*eps0); mu=me*mp_/(me+mp_); a=hbar**2/(mu*k)
def sigb(K): Kj=K*e; rK=(hbar**2*k/(mu*Kj*Kj))**(1/3); return 0.5*math.sqrt(rK*(rK+k/Kj))
def kram(K):
    eta2=B/K; ns=np.arange(1,20000); return ((32*np.pi/(3*np.sqrt(3)))*alpha**3*a*a*eta2**2/(ns*(eta2+ns*ns))).sum()
pool=[]; w=[]
for l in open(err):
    if not l.startswith('FORCEDEVENT'): continue
    d=dict(x.split('=',1) for x in l.split()[1:])
    if d['captured'] not in ('Para-Positronium','Ortho-Positronium'): continue
    E=float(d['E']); L=float(d['L']); K=float(d['K']); b=float(d['b']); P=float(d['P'])
    if not E<0: continue
    n=max(1,round(math.sqrt(B/abs(E))))
    if n>8: continue
    pool.append((n,L)); w.append(P*2*math.pi*sigb(K)**2*math.exp(b*b/(2*sigb(K)**2))/kram(K))
random.seed(557); picks=random.choices(range(len(pool)),weights=w,k=M)
with open(f'samples-{label}.txt','w') as f:
    for i,j in enumerate(picks): f.write(f"{i} {pool[j][0]} {pool[j][1]:.6f}\n")
from collections import Counter
print(label, 'pool', len(pool), 'n counts in sample', sorted(Counter(pool[j][0] for j in picks).items()))
