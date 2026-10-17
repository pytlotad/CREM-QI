# Audit 500: dwell time at n = 2 (Langer 2p start) before the 2 -> 1 photon, from the emission ledger.  python3 ana500.py
import math
hbar=1.054571817e-34; eV=1.602176634e-19
for name,f,B,tqm in (('H','H.err',13.5983,1.596e-9),('Ps','Ps.err',6.80285,3.19e-9)):
    runs=[]; cur=None
    for l in open(f):
        if l.startswith('START'): cur=[]; runs.append(cur)
        elif l.startswith('CREM_LEDGER ') and cur is not None:
            t=l.split(); bars=[k for k,x in enumerate(t) if x=='|']
            def get(side,key):
                s=t[bars[0 if side=='B' else 1]+2:]; return float(s[s.index(key)+1])
            ei=t.index('singlet') if 'singlet' in t else None
            m1=[float(x) for x in t[bars[0]+2:][t[bars[0]+2:].index('m1')+1:][:3]]
            m2=[float(x) for x in t[bars[0]+2:][t[bars[0]+2:].index('m2')+1:][:3]]
            cur.append((float(t[2]),math.sqrt(B*eV/abs(get('B','E'))),math.sqrt(B*eV/abs(get('A','E'))),get('B','L')/hbar,get('A','L')/hbar,m1,m2))
    dw=[]
    for r in runs:
        if r and abs(r[0][1]-2)<1e-3 and abs(r[0][2]-1)<1e-3:
            m1,m2=r[0][5],r[0][6]; d=sum(a*b for a,b in zip(m1,m2))/math.sqrt(sum(a*a for a in m1)*sum(b*b for b in m2))
            dw.append((r[0][0],r[0][3],r[0][4],d))
    ts=sorted(d[0] for d in dw); n=len(ts); mean=sum(ts)/n
    sd=math.sqrt(sum((x-mean)**2 for x in ts)/(n-1)); at_min=sum(1 for x in ts if x-ts[0]<1e-12)
    print(f'{name}: runs {len(runs)}, first photon 2 -> 1 in {n}; dwell mean {mean*1e9:.3f} ns, sd {sd*1e9:.3f} (not exponential: {at_min}/{n} at the minimum), median {ts[n//2]*1e9:.3f}, range {ts[0]*1e9:.3f}-{ts[-1]*1e9:.3f} ns; QM 2p lifetime {tqm*1e9:.3f} ns; ratio {mean/tqm:.3f}')
    print(f'   L before {min(d[1] for d in dw):.4f}-{max(d[1] for d in dw):.4f}, after {min(d[2] for d in dw):.4f}-{max(d[2] for d in dw):.4f}; moment cosine at the photon: '+' '.join(f'{d[3]:+.2f}' for d in dw))
