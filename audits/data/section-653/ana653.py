# Audit 653: energy at formation (CREM_FORMATION_ENERGY=spins|orbit) against 652.
# usage: python3 ana653.py <653 runs> <652 runs>
import re,sys,os,math
d,base=sys.argv[1],sys.argv[2]; eV=1.602176634e-19; u=1e6/eV
F=r'(-?[0-9.e+-]+)'
def parse(f):
    if not os.path.exists(f): return None
    t=open(f).read(); r={}
    m=re.search(r'BORN Q form U '+F+' Z '+F+' E2P '+F+' Q '+F,t); r['Qf']=float(m.group(4))*u if m else None
    m=re.search(r'BORN Q photon .* Qa '+F+r' .* m (-?\d+) dorb '+F,t)
    r['Qa'],r['m'],r['dl']=(float(m.group(1))*u,int(m.group(2)),float(m.group(3))*u) if m else (None,None,None)
    m=re.search(r'eigen photon k=(\d+) w_k=\S+ w_q=(\S+) jump '+F+' classical '+F,t)
    r['k'],r['wq'],r['N']=(int(m.group(1)),m.group(2),(float(m.group(3))-float(m.group(4)))*u) if m else (None,None,None)
    m=re.search(r'Eg (\S+)',t); r['E']=float(m.group(1))/eV if m else None
    m=re.search(r'formation energy spins (hit|miss) t (\S+)',t); r['hit']=m.group(1)=='hit' if m else None
    m=re.search(r'formation energy orbit shift '+F,t); r['shift']=float(m.group(1))*u if m else 0.0
    m=re.search(r'CREM_LEDGER t (\S+)',t); r['t']=float(m.group(1)) if m else None
    r['img']=(re.findall(r'\(1S\) q=\S+ m=(-?\d+) P_S=(\S+)',t) or [('?','?')])[-1]
    r['ok']=None not in (r['Qf'],r['Qa'],r['N'])
    if r['ok']: r['D']=(r['Qf']+r['shift'])-r['dl']-r['Qa']; r['M']=r['N']-r['dl']
    return r
rms=lambda x: math.sqrt(sum(v*v for v in x)/len(x)) if x else float('nan')
for B in ('1e-6','0.1','1'):
    for c,name in (('O1','o-Ps'),('P1','p-Ps')):
        R0=[parse(f'{base}/B{B}_{c}_L2_{i}.err') for i in range(16)]
        print(f'=== B {B} T, {name}')
        D0=[r['D'] for r in R0 if r and r['ok']]
        print(f'  652      : rms D {rms(D0):8.3f}  rms N {rms([r["N"] for r in R0 if r and r["ok"]]):8.3f}  rms M {rms([r["M"] for r in R0 if r and r["ok"]]):8.3f}  rms Q_form {rms([r["Qf"] for r in R0 if r and r["ok"]]):8.3f} ueV')
        for V in ('spins','orbit'):
            R=[parse(f'{d}/{V}_B{B}_{c}_L2_{i}.err') for i in range(16)]
            pairs=[(a,b) for a,b in zip(R0,R) if a and b and a['ok'] and b['ok']]
            use=[(a,b) for a,b in pairs if V=='orbit' or b['hit']]
            hits=sum(1 for r in R if r and r['hit']) if V=='spins' else None
            D=[b['D'] for a,b in use]; D0u=[a['D'] for a,b in use]
            same=[(a,b) for a,b in pairs if a['k']==b['k'] and a['wq']==b['wq'] and a['img']==b['img']]
            dE=max([abs(a['E']-b['E']) for a,b in same],default=0.0)
            dt=[(b['t']-a['t'])/a['t'] for a,b in pairs if a['t'] and b['t']]
            extra=f'hit {hits}/16, ' if V=='spins' else ''
            ratio=rms(D0u)/rms(D) if rms(D)>0 else float('inf')
            print(f'  {V:6s}   : {extra}photons {len(pairs)}/16; rms D {rms(D):8.3f} (652 on same seeds {rms(D0u):8.3f}, ratio {ratio:6.2f})  rms N {rms([b["N"] for a,b in use]):8.3f}  rms M {rms([b["M"] for a,b in use]):8.3f}'
                  f'  rms Q_form,eff {rms([b["Qf"]+b["shift"] for a,b in use]):.2e}  rms delta {rms([b["dl"] for a,b in use]):7.3f}  rms Q_after {rms([b["Qa"] for a,b in use]):7.3f}')
            print(f'             labels (k, q, image) = 652: {len(same)}/{len(pairs)}, max |line - 652| at equal labels {dE:.2e} eV; time to photon rel. change rms {rms(dt):.2e}')
