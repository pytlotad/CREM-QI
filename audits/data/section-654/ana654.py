# Audit 654: named state energy R (CREM_STATE_ENERGY prints) against the 652 prints.
# usage: python3 ana654.py <654 runs> <652 runs>
import re,sys,os,math
d,base=sys.argv[1],sys.argv[2]; eV=1.602176634e-19; u=1e6/eV
F=r'(-?[0-9.e+-]+)'
def kv(line):
    t=line.split(); return {t[i]:float(t[i+1]) for i in range(2,len(t)-1,2)}
def parse(f):
    t=open(f).read(); r={}
    a=[l for l in t.splitlines() if l.startswith('CREM_STATE_ENERGY form')]
    b=[l for l in t.splitlines() if l.startswith('CREM_STATE_ENERGY photon')]
    r['form']=kv(a[0]) if a else None; r['ph']=kv(b[0]) if b else None
    m=re.search(r'BORN Q form U '+F+' Z '+F+' E2P '+F+' Q '+F,t); r['Qf']=float(m.group(4)) if m else None; r['Zf']=float(m.group(2)) if m else None
    m=re.search(r'BORN Q photon .* Zb '+F+' Za '+F+r' .* Qa '+F+r' .* dorb '+F,t)
    r['Zb'],r['Za'],r['Qa'],r['dl']=map(float,m.groups()) if m else (None,)*4
    m=re.search(r'eigen photon k=\d+ w_k=\S+ w_q=\S+ jump '+F+' classical '+F,t); r['N']=float(m.group(1))-float(m.group(2)) if m else None
    return r
rms=lambda x: math.sqrt(sum(v*v for v in x)/len(x)) if x else float('nan')
mx=lambda x: max((abs(v) for v in x),default=float('nan'))
for B in ('1e-6','0.1','1'):
    for c,name in (('O1','o-Ps'),('P1','p-Ps')):
        print(f'=== B {B} T, {name}')
        rows=[]; same=0; n=0
        for i in range(16):
            f=f'{d}/B{B}_{c}_L2_{i}.err'; g=f'{base}/B{B}_{c}_L2_{i}.err'
            a=[l for l in open(f).read().splitlines() if not l.startswith('CREM_STATE_ENERGY')]
            n+=1; same+=a==open(g).read().splitlines()
            r=parse(f); q=parse(g)
            if r['form'] and r['ph'] and q['Qf'] is not None and q['N'] is not None: rows.append((r,q))
        print(f'  H654d other lines = 652: {same}/{n};  photons with R prints {len(rows)}/16')
        # H654a: engine R_after vs offline from 652 prints; named closure
        off=[(-q['Qf']-(q['Zb']-q['Zf'])-q['N']+q['Zb']-q['Za']) for r,q in rows]
        ea=[(r['ph']['Ra']-o)*u for (r,q),o in zip(rows,off)]
        print(f'  H654a R_after engine - offline(652): max {mx(ea):.2e} ueV;  R_form + Q_form max {mx([(r["form"]["R"]+q["Qf"])*u for r,q in rows]):.2e} ueV')
        # D with R named: (Q_form + R_form) - delta - (Q_after + R_after) = (R_def - R_after) - delta  (Q_after = -R_def)
        Dn=[((q['Qf']+r['form']['R'])-q['dl']-(q['Qa']+r['ph']['Ra']))*u for r,q in rows]
        D0=[(q['Qf']-q['dl']-q['Qa'])*u for r,q in rows]
        dWn=[(r['ph']['W']-r['form']['W'])*u for r,q in rows]
        print(f'  D (652) rms {rms(D0):.3f}  ->  D with R rms {rms(Dn):.3f}, max {mx(Dn):.3f} ueV;  D with R + wave work dW: max {mx([a+b for a,b in zip(Dn,dWn)]):.3f} ueV')
        # H654b: carried R_after vs state definition
        db=[(r['ph']['Ra']-r['ph']['Rdef'])*u for r,q in rows]
        print(f'  H654b R_after - R_def: rms {rms(db):.3f}  max {mx(db):.3f} ueV;  rms R_before {rms([r["ph"]["Rb"]*u for r,q in rows]):.3f}, rms R_after {rms([r["ph"]["Ra"]*u for r,q in rows]):.3f}, rms R_def {rms([r["ph"]["Rdef"]*u for r,q in rows]):.3f}')
        # H654c: transport conservation of E_orb + U - W - pending
        X=lambda s,U: s['Eorb']+s[U]-s['W']-s['pending']
        tc=[(X(r['ph'],'U')-X(r['form'],'U'))*u for r,q in rows]
        dU=[(r['ph']['U']-r['form']['U'])*u for r,q in rows]
        gap=[(r['ph']['U']-r['ph']['Uprev'])*u for r,q in rows]
        dW=[(r['ph']['W']-r['form']['W'])*u for r,q in rows]
        print(f'  H654c transport E_orb mu + U - W - pending (photon - formation): rms {rms(tc):.3f} max {mx(tc):.3f} ueV')
        print(f'        dU_cl (orbit) rms {rms(dU):.3f};  U before photon - previousDipoleEnergy rms {rms(gap):.2e};  wave work dW rms {rms(dW):.3f} ueV')
        print(f'        H654b - dU_cl: max |(R_after - R_def) - dU| {mx([a-b for a,b in zip(db,dU)]):.3f} ueV')
        # full ledger (post hoc): D with R + dW against the transport residual, per seed
        print(f'        post hoc: max |(D with R + dW) - (- transport residual)| {mx([a+b+c for a,b,c in zip(Dn,dWn,tc)]):.2e}, max |(D with R + dW) - transport| {mx([a+b-c for a,b,c in zip(Dn,dWn,tc)]):.2e} ueV')
