# Audit 666: the 661 projection import split as B + S2P* - S1S* + dZs.  usage: python3 ana666.py <dir with runs/ b0/> <662 dir>
import re,sys,math,os
u=1.602176634e-25          # J per ueV
d,old=sys.argv[1],sys.argv[2]
K=('Hbar','Hu','before','after','B','S2P','S1S','Zsb','Zsa')
def rd(f):
    if not os.path.exists(f): return None
    t=open(f).read()
    m=re.search(r'PAIR eigen ledger '+' '.join(k+r' (\S+)' for k in K),t)
    i=re.search(r'PAIR eigen photon .* import (\S+)',t)
    if not (m and i): return None
    r={k:float(v)/u for k,v in zip(K,m.groups())}; r['imp']=float(i.group(1))/u
    r['S2Ps']=r['S2P']-r['Zsb']; r['S1Ss']=r['S1S']-r['Zsa']; r['dZs']=r['Zsb']-r['Zsa']
    return r
st=lambda v: (sum(v)/len(v), math.sqrt(sum(x*x for x in v)/len(v)), max(abs(x) for x in v))
allB=[]; okA=okB=True
for sub,ph,tag in (('runs',1,'1.0329 T PHEN 1'),('runs',2,'1.0329 T PHEN 2'),('b0',1,'1e-6 T PHEN 1')):
    R=[r for r in (rd(f'{d}/{sub}/A{ph}_{i}.err') for i in range(64)) if r]
    if not R: print(f'=== {tag}: none'); continue
    idok=max(abs(r['imp']-(r['B']+r['S2Ps']-r['S1Ss']+r['dZs'])) for r in R)
    print(f'=== {tag}: {len(R)} photons; identity import = B + S2P* - S1S* + dZs to {idok:.2e} ueV')
    for k in ('imp','B','S2Ps','S1Ss','dZs'):
        m,q,x=st([r[k] for r in R]); print(f'    {k:5s} mean {m:+9.3f}  rms {q:8.3f}  max|.| {x:8.3f} ueV')
    m,q,x=st([r['Hbar']-r['Hu'] for r in R]); print(f'    Hbar-Hu mean {m:+.3f} rms {q:.3f} max {x:.3f} ueV')
    n2=sum(1 for r in R if abs(r['S2Ps'])<=1); n1=sum(1 for r in R if abs(r['S1Ss'])<=1)
    n0=sum(1 for r in R if abs(r['imp'])<=1e-3)
    print(f'    |S2P*| <= 1 ueV: {n2}/{len(R)};  |S1S*| <= 1 ueV: {n1}/{len(R)};  |import| <= 1e-3 ueV (classical state on the drawn level): {n0}/{len(R)}')
    print(f'    distinct S1S* values (0.01 ueV): {sorted(set(round(r["S1Ss"],2) for r in R))[:8]}')
    okA&=n2==len(R); okB&=n1==len(R); allB+=[r['B'] for r in R]
m=sum(allB)/len(allB); s=math.sqrt(sum((x-m)**2 for x in allB)/(len(allB)-1)/len(allB))
print(f'=== H666a (|S2P*| <= 1 ueV in every photon): {okA};  H666b (|S1S*| <= 1 ueV): {okB};  H666c: mean B {m:+.3f} +- {s:.3f} ueV over {len(allB)}: {abs(m)<=2*s}')
strip=lambda b: b''.join(l for l in b.splitlines(True) if not l.startswith(b'PAIR eigen ledger'))
for sub,ph in (('runs',1),('runs',2),('b0',1)):
    n=sum(1 for i in range(64) if os.path.exists(f'{d}/{sub}/A{ph}_{i}.err') and
          strip(open(f'{d}/{sub}/A{ph}_{i}.err','rb').read())==open(f'{old}/{sub}/A{ph}_{i}.err','rb').read())
    print(f'  H666d {sub}/A{ph}: = 662 byte for byte without the new ledger line {n}/64')
