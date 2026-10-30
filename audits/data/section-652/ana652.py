# Audit 652: where the energy import of 651 comes from (diagnostic prints, engine unchanged).
import re,sys,os,statistics as st
d=sys.argv[1]; ref=sys.argv[2] if len(sys.argv)>2 else None; eV=1.602176634e-19; u=1e6/eV   # J -> ueV
F=r'(-?[0-9.e+-]+)'
def parse(f):
    t=open(f).read(); r={'f':os.path.basename(f)}
    m=re.search(r'BORN Q form U '+F+' Z '+F+' E2P '+F+' Q '+F,t)
    r['form']=tuple(map(float,m.groups())) if m else None
    m=re.search(r'BORN Q photon Ub '+F+' Ua '+F+' Zb '+F+' Za '+F+' E2P '+F+' E1S '+F+' Qb '+F+' Qa '+F+' before '+F+' after '+F+r' m (-?\d+) dorb '+F,t)
    r['ph']=dict(zip('Ub Ua Zb Za E2P E1S Qb Qa before after m dorb'.split(),[float(x) for x in m.groups()])) if m else None
    m=re.search(r'eigen photon k=(\d+) w_k=\S+ w_q=\S+ jump '+F+' classical '+F+' import '+F+' gapShift '+F,t)
    r['eig']=dict(k=int(m.group(1)),jump=float(m.group(2)),cl=float(m.group(3)),imp=float(m.group(4)),gs=float(m.group(5))) if m else None
    r['img']=(re.findall(r'\(1S\) q=\S+ m=(-?\d+) P_S=(\S+)',t) or [('?','?')])[-1]
    return r
def stat(x): return f'mean|.| {st.mean(abs(v) for v in x):8.3f}  min {min(x):9.3f}  max {max(x):9.3f}' if x else 'n/a'
for B in ('1e-6','0.1','1'):
    for c,name in (('O1','o-Ps'),('P1','p-Ps')):
        R=[parse(f'{d}/B{B}_{c}_L2_{i}.err') for i in range(16)]
        R=[r for r in R if r['ph'] and r['eig']]
        print(f'=== B {B} T, {name} ({c}): photons with eigen line and Q print {len(R)}/16')
        I=[r['eig']['imp']*u for r in R]; N=[(r['eig']['jump']-r['eig']['cl'])*u for r in R]
        dl=[r['ph']['dorb']*u for r in R]
        clo=[]
        for r in R:
            p=r['ph']; Rr=(p['Ub']-p['Ua'])-(p['before']-p['after'])+r['eig']['gs']
            clo.append(abs(r['eig']['imp']-((p['Qa']-p['Qb'])+(p['Zb']-p['Za'])+Rr))*u)
        print(f'  H652a closure |I - [(Qa-Qb)+(Zb-Za)+R]| max {max(clo) if clo else 0:.2e} ueV')
        print(f'  I (651 import, ueV)            {stat(I)}')
        print(f'  gapShift (paid by orbit, ueV)  {stat([r["eig"]["gs"]*u for r in R])}')
        print(f'  N = jump - classical (ueV)     {stat(N)}')
        print(f'  delta = orbit - rung (ueV)     {stat(dl)}')
        print(f'  H652f N - delta (ueV)          {stat([a-b for a,b in zip(N,dl)])}')
        Qf=[r['form'][3]*u for r in R if r['form']]; Zf=[r['form'][1]*u for r in R if r['form']]
        print(f'  H652c Q_form (ueV)             {stat(Qf)}   (form prints {len(Qf)})')
        print(f'        Q_before (ueV)           {stat([r["ph"]["Qb"]*u for r in R])}')
        if len(Qf)==len(R): print(f'        Q_before - Q_form (ueV)  {stat([r["ph"]["Qb"]*u-q for r,q in zip(R,Qf)])}')
        for m in (1,-1,0):
            Qa=[r['ph']['Qa']*u for r in R if int(r['ph']['m'])==m]
            print(f'  H652b Q_after image m={m:+d} ({len(Qa):2d})  {stat(Qa)}')
        print(f'  H652d Z_form (ueV)             {stat(Zf)}')
        if len(Zf)==len(R): print(f'        Z_before - Z_form (ueV)  {stat([r["ph"]["Zb"]*u-z for r,z in zip(R,Zf)])}')
        print(f'        Z_before - Z_after (ueV) {stat([(r["ph"]["Zb"]-r["ph"]["Za"])*u for r in R])}')
        if ref:
            same=0; n=0
            for r in R:
                g=f'{ref}/{r["f"]}'
                if not os.path.exists(g): continue
                a=[l for l in open(f'{d}/{r["f"]}').read().splitlines() if not l.startswith('BORN Q ')]
                b=open(g).read().splitlines(); n+=1; same+=a==b
            print(f'  H652e without Q prints = 651 v2: {same}/{n}')
