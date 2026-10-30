# Audit 655: level energy of the 2P label (U_k) against 654.  usage: python3 ana655.py <655 runs> <654 runs>
import re,sys,math
d,base=sys.argv[1],sys.argv[2]; eV=1.602176634e-19; u=1e6/eV
kv=lambda l:(lambda t:{t[i]:float(t[i+1]) for i in range(2,len(t)-1,2)})(l.split())
def parse(f):
    t=open(f).read(); L=t.splitlines(); r={}
    a=[l for l in L if l.startswith('CREM_STATE_ENERGY form')]; b=[l for l in L if l.startswith('CREM_STATE_ENERGY photon')]
    r['form']=kv(a[0]) if a else None; r['ph']=kv(b[0]) if b else None
    m=re.search(r'BORN Q photon .* E2P (\S+) E1S .* dorb (\S+)',t); r['E2Pph'],r['dl']=(float(m.group(1)),float(m.group(2))) if m else (None,None)
    m=re.search(r'eigen photon k=(\d+) w_k=\S+ w_q=(\S+)',t); r['k'],r['wq']=(m.group(1),m.group(2)) if m else (None,None)
    r['img']=(re.findall(r'\(1S\) q=\S+ m=(-?\d+) P_S=(\S+)',t) or [('?','?')])[-1]
    m=re.search(r'CREM_LEDGER t \S+ Eg (\S+)',t); r['Eg']=float(m.group(1)) if m else None
    m=re.search(r'tau (\S+)',t); r['tau']=float.fromhex(m.group(1)) if m else None
    m=re.search(r'CREM_LEDGER t (\S+)',t); r['t']=float(m.group(1)) if m else None
    r['rest']=[l for l in L if not (l.startswith('CREM_STATE_ENERGY') or l.startswith('BORN eigen photon') or l.startswith('BORN Q '))]
    return r
rms=lambda x: math.sqrt(sum(v*v for v in x)/len(x)) if x else float('nan')
mx=lambda x: max((abs(v) for v in x),default=float('nan'))
for B in ('1e-6','0.1','1'):
    for c,name in (('O1','o-Ps'),('P1','p-Ps')):
        R=[parse(f'{d}/B{B}_{c}_L2_{i}.err') for i in range(16)]; R0=[parse(f'{base}/B{B}_{c}_L2_{i}.err') for i in range(16)]
        ok=[(a,b) for a,b in zip(R,R0) if a['ph'] and b['ph']]
        print(f'=== B {B} T, {name}: photons {len(ok)}/16')
        Dn=[((a['ph']['Rdef']-a['ph']['Ra'])-a['dl'])*u for a,b in ok]
        D4=[((b['ph']['Rdef']-b['ph']['Ra'])-b['dl'])*u for a,b in ok]
        print(f'  H655a D with R: rms {rms(Dn):.3e} max {mx(Dn):.3e} ueV   (654: rms {rms(D4):.3f} max {mx(D4):.3f})')
        print(f'        R_after - R_def max {mx([(a["ph"]["Ra"]-a["ph"]["Rdef"])*u for a,b in ok]):.2e};  delta max {mx([a["dl"]*u for a,b in ok]):.2e};'
              f'  U - Uprev max {mx([(a["ph"]["U"]-a["ph"]["Uprev"])*u for a,b in ok]):.2e};  wave work max {mx([(a["ph"]["W"]-a["form"]["W"])*u for a,b in ok]):.2e};'
              f'  E_orb(photon) - E_orb(form) max {mx([(a["ph"]["Eorb"]-a["form"]["Eorb"])*u for a,b in ok]):.2e} ueV')
        print(f'        U_k - (E2P at photon - <Z>_k) max {mx([(a["form"]["Uk"]-(a["E2Pph"]-a["form"]["Zk"]))*u for a,b in ok]):.2e} ueV;  rms U_k - U_cl(form) {rms([(a["form"]["Uk"]-a["form"]["U"])*u for a,b in ok]):.3f} ueV')
        sk=sum(a['k']==b['k'] for a,b in ok); same=[(a,b) for a,b in ok if a['k']==b['k'] and a['wq']==b['wq'] and a['img']==b['img']]
        dE=mx([a['Eg']/b['Eg']-1 for a,b in same])
        print(f'  H655b k = 654: {sk}/{len(ok)};  (k, q, image) = 654: {len(same)}/{len(ok)};  max |line/654 - 1| at equal labels {dE:.2e}')
        dt=[a['t']/b['t']-1 for a,b in ok]; dtau=[a['tau']/b['tau']-1 for a,b in ok]
        print(f'  H655c time to photon rel. change rms {rms(dt):.2e} max {mx(dt):.2e};  tau rel. change rms {rms(dtau):.2e};  other lines = 654: {sum(a["rest"]==b["rest"] for a,b in ok)}/{len(ok)}')
