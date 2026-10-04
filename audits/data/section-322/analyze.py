import sys, math
B=6.80285*1.602176634e-19
for f,bohr in ((sys.argv[1],True),(sys.argv[2],False)):
    try: rows=[list(map(float,l.split()[1:])) for l in open(f) if l.startswith("CREM_REACH") and len(l.split())==15]
    except FileNotFoundError: continue
    first=[r for r in rows if r[11]>1.99 and abs(r[12]-1)<5e-3]
    print(f"{f.split('/')[-1]}: REACH {len(rows)}, first photons {len(first)}")
    if not first: continue
    at1=below=trim=0; worst=0.0; byk={}
    for r in first:
        k=3 if bohr else int(r[10]); Lp=r[12]*r[1]; n=r[11]/math.sqrt(r[0])
        n0=r[11]
        gap=(1-1/n0**2) if bohr else k*2/n0**3          # photon demand in units of B
        nfree=1/math.sqrt(1/n0**2+gap)
        pred=max(1.0 if bohr else nfree, Lp)              # floor only in the Bohr run
        worst=max(worst,abs(n-pred))
        at1+=abs(n-1)<1e-3; below+=n<1-1e-3
        d=byk.setdefault(int(r[10]),[0,0,0]); d[0]+=1; d[1]+=abs(n-1)<1e-3; d[2]+=n<1-1e-3
    ns=sorted(r[11]/math.sqrt(r[0]) for r in first)
    print(f"  n_E after: min {ns[0]:.4f} median {ns[len(ns)//2]:.4f} max {ns[-1]:.4f}")
    print(f"  at n=1 (+-1e-3): {at1}/{len(first)}; below 1: {below}; max |n - max(2/sqrt(1+k), L')| = {worst:.2e}")
    print(f"  by harmonic k: "+", ".join(f"k={k}: {v[0]} (n=1: {v[1]}, <1: {v[2]})" for k,v in sorted(byk.items())))
    print(f"  L' quartiles: "+" ".join(f"{x:.3f}" for x in [sorted(r[12]*r[1] for r in first)[i*len(first)//4] for i in (1,2,3)]))
