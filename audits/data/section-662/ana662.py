# Audit 662: prepared singlet projected as |L^> (x) |S> at the photon (CREM_PAIR_FIELD_EIGEN=1).  usage: python3 ana662.py <dir with b0/ runs/> <661 dir>
import re,sys,math,os
h=6.62607015e-34; c0=299792458.0; e0=1.602176634e-19
line=lambda f: f'{f:.1f} MHz ({c0/(f*1e6)*1e9:.7f} nm; {h*f*1e6/e0:.9f} eV)'
d,old=sys.argv[1],sys.argv[2]
def img(f):
    if not os.path.exists(f): return None
    t=open(f).read()
    m=re.search(r'PAIR field eigen \(1S\) q=\S+ m=(-?\d+) P_S=(\S+)',t); g=re.search(r'CREM_LEDGER t \S+ Eg (\S+)',t)
    if not (m and g): return None
    return dict(m=int(m.group(1)),PS=float(m.group(2)),E=float(g.group(1))/h/1e6)
def binom(k,n,p): s=math.sqrt(p*(1-p)/n); return k/n,s,abs(k/n-p)<=2*s
for tag,dd,p0 in (('662 (singlet |S>)',d,None),('661 (product state)',old,None)):
    R=[r for r in (img(f'{dd}/b0/A1_{i}.err') for i in range(64)) if r]
    if not R: print(f'=== {tag}: no zero-field PHEN 1 probes'); continue
    k=sum(1 for r in R if r['PS']>=0.5); f,s,ok=binom(k,len(R),5/9)
    print(f'=== {tag}, BEXT=1e-6 PHEN 1: {len(R)} photons, 1S F=0 {k}/{len(R)} = {f:.3f} (rule = QM 0.556 +- {2*s:.3f} at 2 sigma): H662b {ok}')
    for v in sorted(set(round(r['E'],1) for r in R)):
        n=sum(1 for r in R if round(r['E'],1)==v); F=[r['PS']>=0.5 for r in R if round(r['E'],1)==v]
        print(f'    {line(v)}: {n} photons, 1S F={0 if F[0] else 1}')
R=[r for r in (img(f'{d}/runs/A1_{i}.err') for i in range(64)) if r]
k=sum(1 for r in R if r['m']==-1); f,s,ok=binom(k,len(R),0.0602)
ku=sum(1 for r in R if r['m']==1)
print(f'=== 662 BEXT=1.0329 PHEN 1: {len(R)} photons; 1S |down down> (1Sd) {k}/{len(R)} = {f:.3f} (rule = QM 0.060 +- {2*s:.3f}): H662c {ok};'
      f' |up up> {ku}/{len(R)} = {ku/len(R):.3f} (0.060); m = 0 {len(R)-k-ku}/{len(R)} (0.880)')
R1=[r for r in (img(f'{old}/runs/A1_{i}.err') for i in range(64)) if r]
lines1=set(round(r['E'],1) for r in R1); lines2=set(round(r['E'],1) for r in R)
print(f'  lines of 662 PHEN 1 that 661 (both PHEN, 128 probes) also has: {len(lines2 & set(round(r["E"],1) for ph in (1,2) for r in (img(f"{old}/runs/A{ph}_{i}.err") for i in range(64)) if r))}/{len(lines2)}')
# the only new stderr line of 662 is the CREM_DEBUG_BORN print 'PAIR prepared singlet a -> b'
strip=lambda b: b''.join(l for l in b.splitlines(True) if not l.startswith(b'PAIR prepared singlet'))
same=sum(1 for i in range(64) if os.path.exists(f'{d}/runs/A2_{i}.err') and strip(open(f'{d}/runs/A2_{i}.err','rb').read())==open(f'{old}/runs/A2_{i}.err','rb').read())
print(f'  PHEN 2 at 1.0329 T = 661 byte for byte (without the new PAIR prepared singlet line): {same}/64')
