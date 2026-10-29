# Audit 649: field labels at formation, off-rung starts vs the Langer start (per cell B, config, seed).
import re,sys,glob,os
d=sys.argv[1]
def parse(f):
    t=open(f).read()
    born=[l for l in t.splitlines() if l.startswith('BORN')]
    eg=[float(x) for x in re.findall(r'CREM_LEDGER t \S+ Eg (\S+)',t)]
    m=re.search(r'stop (\d+) photons (\d+) tau (\S+) born (\S+) chosen (\d)',t)
    return born,eg,float.fromhex(m.group(3)),int(m.group(5))
pairs=[('c1','L1'),('c2','L2'),('n23','L2'),('prec1','L1'),('prec2','L2'),('pren23','L2')]
stat={}
for B in ['0.01','0.1','1']:
  for c in ['Q0','P1','O1']:
    for i in range(4):
      ref={k:parse(f'{d}/B{B}_{c}_{k}_{i}.err') for k in ['L1','L2']}
      for a,b in pairs:
        A=parse(f'{d}/B{B}_{c}_{a}_{i}.err'); R=ref[b]
        same_born=A[0]==R[0]; same_n=len(A[1])==len(R[1])
        de=max([abs(x/y-1) for x,y in zip(A[1],R[1])],default=0.0) if same_n else float('inf')
        dt=A[2]/R[2]-1
        ch='p-Ps' if R[3]==1 or c=='P1' else 'o-Ps'
        if c=='Q0': ch=('p-Ps' if R[3] else 'o-Ps')+'(Born)'
        key=(a,c); s=stat.setdefault(key,[0,0,0.0,0.0,[]]); s[0]+=1; s[1]+=same_born and same_n and de<1e-12
        s[2]=max(s[2],de if de!=float('inf') else 9.9); s[3]=max(s[3],abs(dt))
        if not same_born or not same_n or de>=1e-12: s[4].append(f'B{B} i{i} {ch}: '+' / '.join(A[0])+'  vs  '+' / '.join(R[0])+f'  dE {de:.2e}')
        print(f'B{B} {c} i{i} {a:6s} vs {b} [{ch}] born_same {int(same_born)} nphot {len(A[1])}/{len(R[1])} dE {de:.2e} dtau {dt:+.3e}')
print('\nSUMMARY (start, config): cells, identical (Born lines + photon energies 1e-12), max dE, max |dtau|')
for k,s in stat.items():
    print(f'{k[0]:6s} {k[1]}: {s[1]}/{s[0]}  dE {s[2]:.2e}  |dtau| {s[3]:.3e}')
    for x in s[4][:6]: print('     ',x)
