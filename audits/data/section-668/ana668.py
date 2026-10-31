# Audit 668: projection inside the classical j (CREM_PAIR_J_PROJECTION=1).  usage: python3 ana668.py <dir with b0/ runs/>
import re,sys,math,os
h=6.62607015e-34; c0=299792458.0; e0=1.602176634e-19; u=1.602176634e-25
line=lambda f: f'{f:.1f} MHz ({c0/(f*1e6)*1e9:.7f} nm; {h*f*1e6/e0:.9f} eV)'
d=sys.argv[1]
def rd(f):
    if not os.path.exists(f): return None
    t=open(f).read()
    if not re.search(r'^END',t,re.M): return None
    i=re.search(r'PAIR eigen photon .* import (\S+)',t); g=re.search(r'CREM_LEDGER t \S+ Eg (\S+)',t)
    q=re.search(r'PAIR field eigen \(1S\) q=\S+ m=(-?\d+) P_S=(\S+)',t)
    L=re.search(r'PAIR eigen ledger .* S2P (\S+) S1S (\S+) Zsb (\S+) Zsa (\S+)',t); j=re.search(r'PAIR j projection jcl (\S+) kept (\S+)',t)
    if not (i and g and q and L): return None
    return dict(imp=float(i.group(1))/u,E=float(g.group(1))/h/1e6,m=int(q.group(1)),PS=float(q.group(2)),
                S2P=(float(L.group(1))-float(L.group(3)))/u,S1S=(float(L.group(2))-float(L.group(4)))/u,
                jcl=float(j.group(1)) if j else None,kept=float(j.group(2)) if j else None)
def binom(k,n,p): s=math.sqrt(p*(1-p)/max(n,1)); return k/max(n,1),s,abs(k/max(n,1)-p)<=2*s
st=lambda v: (max(abs(x) for x in v), math.sqrt(sum(x*x for x in v)/len(v)))
R=[r for r in (rd(f'{d}/b0/A1_{i}.err') for i in range(64)) if r]
if R:
    mx,rms=st([r['imp'] for r in R]); k=sum(1 for r in R if r['PS']>=0.5); f,s,ok=binom(k,len(R),5/9)
    n32=sum(1 for r in R if r['jcl']==1.5)
    print(f'=== 1e-6 T PHEN 1: {len(R)} photons; |import| max {mx:.4f} ueV, rms {rms:.4f}; <= 1 ueV in {sum(1 for r in R if abs(r["imp"])<=1)}/{len(R)}')
    print(f'    1S F=0 {k}/{len(R)} = {f:.3f} (QM 0.556 +- {2*s:.3f}; 667 prediction 0.531): H668a {ok and mx<=1}; classical j = 3/2 in {n32}/{len(R)}')
    for v in sorted(set(round(r['E'],1) for r in R)):
        n=sum(1 for r in R if round(r['E'],1)==v); print(f'    {line(v)}: {n} photons')
for ph in (1,2):
    R=[r for r in (rd(f'{d}/runs/A{ph}_{i}.err') for i in range(64)) if r]
    if not R: continue
    mx,rms=st([r['S2P'] for r in R]); mi,ri=st([r['imp'] for r in R])
    k=sum(1 for r in R if r['m']==-1)
    print(f'=== 1.0329 T PHEN {ph}: {len(R)} photons; |S2P*| max {mx:.3f} rms {rms:.3f} ueV, <= 1 ueV in {sum(1 for r in R if abs(r["S2P"])<=1)}/{len(R)};'
          f' |import| max {mi:.3f} rms {ri:.3f} ueV; 1Sd {k}/{len(R)} = {k/len(R):.3f}; classical j = 3/2 in {sum(1 for r in R if r["jcl"]==1.5)}/{len(R)}')
