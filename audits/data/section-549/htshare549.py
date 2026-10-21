# Audit 549: level shares of captures weighted by Horvitz-Thompson weights of the Rayleigh-sampled b (w ~ 2 pi b / p(b), p(b) = b/s^2 exp(-b^2/2s^2)), so the shares are the
# cross-section shares sigma_n/sigma (raw counts favour small b, hence n = 1).  python3 htshare549.py CAPDUMPfile sigma_b_pm B_eV
import sys, re, math
f,s,B=sys.argv[1],float(sys.argv[2])*1e-12,float(sys.argv[3]); W={}; tot=0
for l in open(f):
    if not l.startswith('CAPDUMP'): continue
    out=re.search(r'out=(\S+)',l).group(1)
    if out in ('Scattering','Collision') or 'nresolved' in out: continue
    E=float(re.search(r'Ecoul=\S+/(\S+)',l).group(1)); b=float(re.search(r' b=(\S+)',l).group(1))
    if E>=0: continue
    n=max(1,round(math.sqrt(B/abs(E)))); w=2*math.pi*s*s*math.exp(b*b/(2*s*s)); W[n]=W.get(n,0)+w; tot+=w
print("HT-weighted shares: "+' '.join(f"n={n}: {W[n]/tot:.3f}" for n in sorted(W)))
