# Audit 650: field line distribution of the Langer start vs path starts (o-Ps O1: L2, c2, n23; p-Ps P1: L2), B = 0.1 and 1 T, seeds 0-31.
import re,sys,math
d=sys.argv[1]; eV=1.602176634e-19; N=32
def parse(f):
    t=open(f).read()
    born=[l for l in t.splitlines() if l.startswith('BORN')]
    eg=[float(x)/eV for x in re.findall(r'CREM_LEDGER t \S+ Eg (\S+)',t)]
    m=re.search(r'stop (\d+) photons (\d+)',t)
    img=[l.split('-> ')[-1] for l in born if '(1S)' in l]
    m1=re.findall(r'\(1S\) q=\S+ m=(-?\d+)',t)
    return dict(born=born,E=eg[0] if eg else float('nan'),stop=int(m.group(1)),img=(img[-1] if img else '-')+(' m='+m1[-1] if m1 else ''))
def fisher(a,b,c,d_):  # two-sided Fisher exact p for [[a,b],[c,d]]
    n=a+b+c+d_; r1=a+b; c1=a+c
    def p(x): return math.comb(r1,x)*math.comb(n-r1,c1-x)/math.comb(n,c1)
    p0=p(a); return min(1.0,sum(p(x) for x in range(max(0,c1-(n-r1)),min(r1,c1)+1) if p(x)<=p0*(1+1e-9)))
for B in ['0.1','1']:
    R={(c,s):[parse(f'{d}/B{B}_{c}_{s}_{i}.err') for i in range(N)] for c,s in [('O1','L2'),('O1','c2'),('O1','n23'),('P1','L2')]}
    allE=[r['E'] for v in R.values() for r in v]
    D=sorted({x for x in allE if sum(abs(y/x-1)<1e-12 for y in allE)>=2})
    dev=lambda E: min(abs(E-x) for x in D)
    print(f'=== B = {B} T: discrete set D ({len(D)} values, eV):',' '.join(f'{x:.10f}' for x in D))
    for k,v in R.items():
        off=[r for r in v if dev(r['E'])/r['E']>1e-12]
        print(f'  {k[0]} {k[1]:3s}: stops {sorted({r["stop"] for r in v})}  off-D {len(off)}/{N}  max dev {max([dev(r["E"]) for r in v]):.3e} eV')
        hist={}
        for r in v: hist.setdefault(r['img'],[]).append(r)
        for im,rs in sorted(hist.items()):
            Es=[r['E'] for r in rs]; o=sum(dev(e)/e>1e-12 for e in Es)
            print(f'      1S {im:22s}: {len(rs):2d} runs, off-D {o}, lines {min(Es):.10f} - {max(Es):.10f} (range {max(Es)-min(Es):.3e} eV)')
    L=R[('O1','L2')]; offL=sum(dev(r['E'])/r['E']>1e-12 for r in L); maxdevL=max(dev(r['E']) for r in L)
    for s in ['c2','n23']:
        P=R[('O1',s)]; offP=sum(dev(r['E'])/r['E']>1e-12 for r in P)
        diff=[(i,P[i],L[i]) for i in range(N) if abs(P[i]['E']/L[i]['E']-1)>=1e-12]
        sameborn=[x for x in diff if x[1]['born']==x[2]['born']]
        maxd=max([abs(a['E']-b['E']) for _,a,b in diff],default=0.0)
        sameimg=sum(a['img']==b['img'] for _,a,b in sameborn)
        print(f'  H650b {s} vs L2: (i) off-D {offP}/{N} vs {offL}/{N}, Fisher p = {fisher(offP,N-offP,offL,N-offL):.3f};'
              f' (ii) differing cells {len(diff)} (labels equal {len(sameborn)}), max |dE| {maxd:.3e} eV vs Langer max dev {maxdevL:.3e} eV;'
              f' (iii) same 1S image {sameimg}/{len(sameborn)}')
        for i,a,b in diff: print(f'      seed {i}: {s} {a["E"]:.10f} [{a["img"]}]  L2 {b["E"]:.10f} [{b["img"]}]  born_same {int(a["born"]==b["born"])}')
