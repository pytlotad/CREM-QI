# Audit 518: per-photon energy closure within each run: E_before - E_after - Eg + (UF_before - UF_after) + Delta m hbar omega_L B
# (orbit energies E from the ledger; ~10 MHz background = recoil / pending spin energy, present at B = 0 too).  python3 closure518.py files...
import math, sys
hbar=1.054571817e-34; muB=9.2740100783e-24; mp=1.67262192369e-27; me=9.1093837015e-31; h=6.62607015e-34; wL=muB*mp/(mp+me)
def V(a): return [float(x) for x in a[:3]]
def dot(a,b): return sum(x*y for x,y in zip(a,b))
for f in sys.argv[1:]:
    idx=None; B=None; worst=0
    for l in open(f):
        if l.startswith('BFIELD'): B=V(l.split()[3:6])
        elif l.startswith('START'): idx=int(l.split()[2])
        elif l.startswith('CREM_LEDGER '):
            t=l.split(); b=[k for k,x in enumerate(t) if x=='|']; sb=t[b[0]+2:]; sa=t[b[1]+2:]; g=lambda s,k: float(s[s.index(k)+1])
            Bm=math.sqrt(dot(B,B)); m=round(dot(V(t[t.index('sphoton')+1:]),[x/Bm for x in B])/hbar) if Bm>0 else 0
            c=(g(sb,'E')-g(sa,'E')-float(t[t.index('Eg')+1])+(g(sb,'UF')-g(sa,'UF'))+m*wL*Bm)/h*1e-6
            worst=max(worst,abs(c)); print(f"{f} idx {idx} L {g(sb,'L')/hbar:.3f} -> {g(sa,'L')/hbar:.3f} dm {m:+d} closure {c:+.3f} MHz")
    print(f"{f}: max |closure| {worst:.3f} MHz")
