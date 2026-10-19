# Audit 518: 3d -> 2p in 1 T vs B = 0 (same seeds, both with the current code): Eg(B) - Eg(0) = Delta m hbar omega_L + [UF(B) - UF(0)],
# UF = level energy beyond Bohr released (before - after; with import 21 it holds the j of 3d and 2p, which the field's spin precession can
# change between the paired runs).  Uses ana504's loader.  python3 ana518.py D1T.err D0.err
import sys, math, runpy
L=runpy.run_path('../section-504/ana504.py',run_name='lib') if False else None
hbar=1.054571817e-34; muB=9.2740100783e-24; mp=1.67262192369e-27; me=9.1093837015e-31; h=6.62607015e-34
def V(a): return [float(x) for x in a[:3]]
def dot(a,b): return sum(x*y for x,y in zip(a,b))
def unit(a): n=math.sqrt(dot(a,a)); return [x/n for x in a]
def load(f):
    out={}; idx=None; B=None
    for l in open(f):
        if l.startswith('BFIELD'): B=[float(x) for x in l.split()[3:6]]
        elif l.startswith('START'): idx=int(l.split()[2])
        elif l.startswith('CREM_LEDGER ') and idx is not None and idx not in out:
            t=l.split(); b=[k for k,x in enumerate(t) if x=='|']; sb=t[b[0]+2:]; sa=t[b[1]+2:]; g=lambda s,k: s[s.index(k)+1:]
            out[idx]=dict(B=B,Eg=float(t[t.index('Eg')+1]),sph=V(t[t.index('sphoton')+1:]),UF=float(g(sb,'UF')[0])-float(g(sa,'UF')[0]),
                          UFB=float(g(sb,'UF')[0]),UFA=float(g(sa,'UF')[0]))
    return out
Z=load(sys.argv[1]); O=load(sys.argv[2]); wL=muB*mp/(mp+me); worst=0
for k in sorted(set(Z)&set(O)):
    z,o=Z[k],O[k]; Bm=math.sqrt(dot(z['B'],z['B'])); m=round(dot(z['sph'],unit(z['B']))/hbar)
    d=(z['Eg']-o['Eg'])/(muB*Bm); dUF=(z['UF']-o['UF'])/(muB*Bm); r=d-m*wL/muB-dUF; worst=max(worst,abs(r))
    print(f"  idx {k}: Delta m {m:+d}; (Eg_B - Eg_0)/(mu_B B) {d:+.4f}; Delta m omega_L {m*wL/muB:+.5f}; [UF_B - UF_0]/(mu_B B) {dUF:+.4f}"
          f" (UF after: {z['UFA']/h*1e-6:+.1f} / {o['UFA']/h*1e-6:+.1f} MHz); residual {r:+.2e}")
print(f"max |residual| = {worst:.2e} mu_B B")
