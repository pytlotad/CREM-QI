# Audit 467: engine precession sectors vs Breit-Pauli spin-orbit (SO) and spin-spin (SS), per particle.
import re, math
def load(f):
    out={}; pair=None; cfg=None
    for l in open(f):
        if l.startswith('# ') and ':' in l and ' g ' in l: pair=l.split(':')[0][2:]; continue
        m=re.match(r'\s+\[(.*)\]',l)
        if m: cfg=m.group(1); continue
        p=l.split()
        if len(p)==4 and pair and cfg and p[0] in ('engine1','engine2','SO1','SO2','SS1','SS2','SO+SS1','SO+SS2'):
            out[(pair,cfg,p[0])]=[float(x) for x in p[1:]]
    return out
F,C,D=load('p467-full.txt'),load('p467-charge.txt'),load('p467-dipole.txt')
def n(v): return math.sqrt(sum(x*x for x in v))
def rel(a,b): return n([x-y for x,y in zip(a,b)])/max(n(b),1e-300)
for key in sorted(set(k[:2] for k in F)):
    for i in ('1','2'):
        e='engine'+i
        rc=rel(C[key+(e,)],C[key+('SO'+i,)]); rd=rel(D[key+(e,)],D[key+('SS'+i,)])
        lin=rel(F[key+(e,)],[a+b for a,b in zip(C[key+(e,)],D[key+(e,)])]); rf=rel(F[key+(e,)],F[key+('SO+SS'+i,)])
        print(f'{key[0]:18s} {key[1]:20s} p{i}: charge/SO {rc:.2e}  dipole/SS {rd:.2e}  full/(charge+dipole) {lin:.2e}  full/(SO+SS) {rf:.2e}  |SS|/|SO| {n(D[key+("SS"+i,)])/n(C[key+("SO"+i,)]):.3f}')
