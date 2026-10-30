# Audit 657: level energy from the transport (CREM_PS_TRANSPORT_LEVEL=1) against 654 / 655.
# usage: python3 ana657.py <657 runs> <654 runs> <655 runs> <657 runs0 (B = 0 without the switch)>
import re,sys,math
d,base,b655,d0=sys.argv[1:5]; eV=1.602176634e-19; u=1e6/eV
kv=lambda l:(lambda t:{t[i]:float(t[i+1]) for i in range(2,len(t)-1,2)})(l.split())
def parse(f):
    t=open(f).read(); L=t.splitlines(); r={'lines':L}
    a=[l for l in L if l.startswith('CREM_STATE_ENERGY form')]; b=[l for l in L if l.startswith('CREM_STATE_ENERGY photon')]
    r['form']=kv(a[0]) if a else None; r['ph']=kv(b[0]) if b else None
    m=re.search(r'BORN Q photon .* dorb (\S+)',t); r['dl']=float(m.group(1)) if m else None
    m=re.search(r'eigen photon k=(\d+) w_k=\S+ w_q=(\S+)',t); r['k'],r['wq']=(m.group(1),m.group(2)) if m else (None,None)
    r['img']=(re.findall(r'\(1S\) q=\S+ m=(-?\d+) P_S=(\S+)',t) or [('?','?')])[-1]
    m=re.search(r'CREM_LEDGER t (\S+) Eg (\S+)',t); r['t'],r['Eg']=(float(m.group(1)),float(m.group(2))) if m else (None,None)
    m=re.search(r'tau (\S+)',t); r['tau']=float.fromhex(m.group(1)) if m else None
    m=re.search(r'CREM_TRANSPORT_LEVEL photon C (\S+) sumdF (\S+)',t); r['C'],r['dF']=(float(m.group(1)),float(m.group(2))) if m else (0.0,0.0)
    return r
rms=lambda x: math.sqrt(sum(v*v for v in x)/len(x)) if x else float('nan')
mx=lambda x: max((abs(v) for v in x),default=float('nan'))
D=lambda r: ((r['ph']['Rdef']-r['ph']['Ra'])-r['dl'])*u
for B in ('1e-6','0.1','1'):
    for c,name in (('O1','o-Ps'),('P1','p-Ps')):
        R=[parse(f'{d}/B{B}_{c}_L2_{i}.err') for i in range(16)]; R0=[parse(f'{base}/B{B}_{c}_L2_{i}.err') for i in range(16)]
        R5=[parse(f'{b655}/B{B}_{c}_L2_{i}.err') for i in range(16)]
        ok=[(a,b,e) for a,b,e in zip(R,R0,R5) if a['ph'] and b['ph'] and e['ph']]
        print(f'=== B {B} T, {name}: photons {len(ok)}/16')
        Dn=[D(a) for a,b,e in ok]
        print(f'  H657a D with R: rms {rms(Dn):.3e} max {mx(Dn):.3e} ueV   (654: rms {rms([D(b) for a,b,e in ok]):.3f}; 655 U_k: rms {rms([D(e) for a,b,e in ok]):.1e})')
        print(f'        |D + sum dF_tr| max {mx([D(a)+a["dF"]*u for a,b,e in ok]):.2e} ueV;  sum dF_tr rms {rms([a["dF"]*u for a,b,e in ok]):.2e} max {mx([a["dF"]*u for a,b,e in ok]):.2e};'
              f'  C rms {rms([a["C"]*u for a,b,e in ok]):.3f} max {mx([a["C"]*u for a,b,e in ok]):.3f} ueV')
        sk=sum(a['k']==b['k'] for a,b,e in ok); same=[(a,b) for a,b,e in ok if a['k']==b['k'] and a['wq']==b['wq'] and a['img']==b['img']]
        print(f'  H657b k = 654: {sk}/{len(ok)};  (k, q, image) = 654: {len(same)}/{len(ok)};  max |line/654 - 1| at equal labels {mx([a["Eg"]/b["Eg"]-1 for a,b in same]):.2e}')
        dE=[(a['ph']['Eorb']-a['form']['Eorb'])*u for a,b,e in ok]; dW=[(a['ph']['W']-a['form']['W'])*u for a,b,e in ok]
        dW0=[(b['ph']['W']-b['form']['W'])*u for a,b,e in ok]
        print(f'  H657c E_orb(photon) - E_orb(form) max {mx(dE):.2e};  wave work max {mx(dW):.2e} ueV  (654: wave work rms {rms(dW0):.3f} max {mx(dW0):.3f})')
        dt=[a['t']/b['t']-1 for a,b,e in ok]
        print(f'        time to photon rel. change vs 654 rms {rms(dt):.2e} max {mx(dt):.2e}')
print('=== B = 0 T (classical photon path): with the switch vs without')
for c,name in (('O1','o-Ps'),('P1','p-Ps')):
    R=[parse(f'{d}/B0_{c}_L2_{i}.err') for i in range(16)]; R0=[parse(f'{d0}/B0_{c}_L2_{i}.err') for i in range(16)]
    ok=[(a,b) for a,b in zip(R,R0) if a['Eg'] and b['Eg']]
    byte=sum([l for l in a['lines'] if not l.startswith('CREM_TRANSPORT_LEVEL')]==b['lines'] for a,b in ok)
    print(f'  {name}: photons {len(ok)}/16; H657d max |line rel| {mx([a["Eg"]/b["Eg"]-1 for a,b in ok]):.2e}, max |t rel| {mx([a["t"]/b["t"]-1 for a,b in ok]):.2e},'
          f' max |tau rel| {mx([a["tau"]/b["tau"]-1 for a,b in ok]):.2e}; byte-identical (without the 657 print) {byte}/{len(ok)};'
          f'  C max {mx([a["C"]*u for a,b in ok]):.2e} ueV, sum dF max {mx([a["dF"]*u for a,b in ok]):.2e} ueV')
# post hoc (after the first output): the orbit off the rung at the photon (delta != 0 in some probes: the last dF_tr booked by the
# exchange not yet moved to the wave) enters D twice -- the closure with the orbit's own change.
print('=== post hoc: D + sum dF_tr + (E_orb(photon) - E_orb(form)), o-Ps')
for B in ('1e-6','0.1','1'):
    R=[parse(f'{d}/B{B}_O1_L2_{i}.err') for i in range(16)]
    x=[D(a)+a['dF']*u+(a['ph']['Eorb']-a['form']['Eorb'])*u for a in R if a['ph']]
    print(f'  B {B} T: max {mx(x):.2e} ueV; probes with orbit off the rung at the photon {sum(1 for a in R if a["ph"] and a["dl"]!=0.0)}/16')
