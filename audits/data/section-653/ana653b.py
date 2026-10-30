# Audit 653, post hoc (after the main analysis): (i) non-trivial hits of A (|Q_form| > 1e-35 J before rotation);
# (ii) variant B: does the orbit keep the shift?  delta_B - delta_652 against the shift -Q_form, per seed.
import re,sys,os,math
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
d,base=sys.argv[1],sys.argv[2]; eV=1.602176634e-19; u=1e6/eV
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'ana653.py')).read().split('rms=lambda')[0].split("d,base=sys.argv[1],sys.argv[2]; eV=1.602176634e-19; u=1e6/eV")[1])
rms=lambda x: math.sqrt(sum(v*v for v in x)/len(x)) if x else float('nan')
for B in ('1e-6','0.1','1'):
    for c,name in (('O1','o-Ps'),('P1','p-Ps')):
        print(f'=== B {B} T, {name}')
        nt=0; nth=0
        for i in range(16):
            t=open(f'{d}/spins_B{B}_{c}_L2_{i}.err').read()
            m=re.search(r'formation energy spins (hit|miss) t \S+ Q0 (\S+)',t)
            if abs(float(m.group(2)))>1e-35: nt+=1; nth+=m.group(1)=='hit'
        print(f'  A: non-trivial formations (|Q_form| > 1e-35 J) {nt}/16, of them hit {nth}')
        rows=[]
        for i in range(16):
            a=parse(f'{base}/B{B}_{c}_L2_{i}.err'); b=parse(f'{d}/orbit_B{B}_{c}_L2_{i}.err')
            if a['ok'] and b['ok'] and abs(b['shift'])>1e-6: rows.append((i,b['shift'],a['dl'],b['dl']))
        kept=[(db-da)/s for i,s,da,db in rows]
        print(f'  B: seeds with shift > 1e-6 ueV {len(rows)}; (delta_B - delta_652)/shift: '+' '.join(f'{k:.3f}' for k in kept))
        print(f'     rms shift {rms([s for i,s,da,db in rows]):.3f} ueV, rms delta_652 {rms([da for i,s,da,db in rows]):.3f}, rms delta_B {rms([db for i,s,da,db in rows]):.3f}')
