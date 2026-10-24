# Audit 614: channel of the first photon per run from CREM_LEDGER lines; binomial tests against QM and the offline engine shares.
import glob, math, sys
hbar=1.054571817e-34; B=2.179872361e-18/(1+1/1836.15267343)
def first_channels(pattern):
    res=[]; nophoton=0; stops={}
    for f in sorted(glob.glob(pattern)):
        cur=None; got=False
        for line in open(f):
            if line.startswith('START idx'): cur=line.split()[2]; got=False
            elif line.startswith('CREM_LEDGER t') and not got:
                t=line.split(); a=t.index('A'); Ea=float(t[a+2]); La=float(t[a+4])
                res.append((round(math.sqrt(B/abs(Ea))),round(La/hbar-0.5))); got=True
            elif line.startswith('END idx'):
                st=line.split()[4]; stops[st]=stops.get(st,0)+1
                if not got: nophoton+=1
    return res,nophoton,stops
def binom_p(k,n,p):   # two-sided exact
    from math import comb
    pk=comb(n,k)*p**k*(1-p)**(n-k)
    return min(1.0,sum(comb(n,i)*p**i*(1-p)**(n-i) for i in range(n+1) if comb(n,i)*p**i*(1-p)**(n-i)<=pk*(1+1e-12)))
for tag,pat,tests in (("3p new","runs/3p_new_*.err",[((2,0),0.118,"QM / offline new"),((2,0),0.087,"offline old")]),
                      ("3p old","runs/3p_old_*.err",[((2,0),0.087,"offline old"),((2,0),0.118,"QM")]),
                      ("4p new","runs/4p_new_*.err",[((2,0),0.122,"offline new"),((2,0),0.119,"QM"),((2,0),0.037,"offline old"),((3,0),0.038,"QM 3s")])):
    ch,nop,stops=first_channels(pat); n=len(ch); cnt={}
    for c in ch: cnt[c]=cnt.get(c,0)+1
    print(f"{tag}: {n} first photons, runs without photon {nop}, stop causes {stops}; channels {dict(sorted(cnt.items()))}")
    for c,p,lab in tests:
        k=cnt.get(c,0); print(f"   {c[0]}{'spd'[c[1]]}: {k}/{n} = {k/n:.3f} +- {math.sqrt(p*(1-p)/n):.3f}  vs {p} ({lab}): p = {binom_p(k,n,p):.3g}")
