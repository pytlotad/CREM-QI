# Audit 516: production check of import 21 (default) -- H 2p -> 1s photons from the emission ledger, default vs CREM_NO_SOMMERFELD_J=1.
# Per photon: t, Eg, the level energy beyond Bohr before (UF_B), j read from the before-state (cos(L, s_e) vs -1/3, s_e = -m2^), and
# Eg(default) - Eg(off) for the same run index.  python3 ana516.py def.err off.err
import sys, math
h=6.62607015e-34; MHz=1e-6/h
def rows(f):
    out=[]
    for line in open(f):
        if not line.startswith('CREM_LEDGER '): continue
        t=line.split(); bars=[k for k,x in enumerate(t) if x=='|']; head=t[:bars[0]]; B=t[bars[0]+1:bars[1]]; A=t[bars[1]+1:]
        g=lambda seg,key,n=1: [float(x) for x in seg[seg.index(key)+1:seg.index(key)+1+n]]
        ax=g(B,'ax',3); m2=g(B,'m2',3)
        c=-sum(a*b for a,b in zip(ax,m2))/math.sqrt(sum(a*a for a in ax)*sum(b*b for b in m2))
        out.append(dict(t=g(head,'t')[0],Eg=g(head,'Eg')[0],UFB=g(B,'UF')[0],UFA=g(A,'UF')[0],LB=g(B,'L')[0]/1.054571817e-34,c=c,EB=g(B,'E')[0],EA=g(A,'E')[0]))
    return out
D=rows(sys.argv[1]); O=rows(sys.argv[2])
print(f"photons default {len(D)}, off {len(O)} (per run: 2p -> 1s, then 1s stop)")
for i,(d,o) in enumerate(zip(D,O)):
    if d['LB']<1: continue
    j='3/2' if d['c']>-1/3 else '1/2'
    print(f"  photon {i}: L_before {d['LB']:.4f} hbar, cos(L,s_e) {d['c']:+.4f} -> j {j}; t {d['t']*1e9:.4f} / {o['t']*1e9:.4f} ns; UF_before {d['UFB']*MHz:+.2f} / {o['UFB']*MHz:+.2f} MHz;"
          f" UF_after {d['UFA']*MHz:+.2f} / {o['UFA']*MHz:+.2f}; Eg(default) - Eg(off) {(d['Eg']-o['Eg'])*MHz:+.2f} MHz;"
          f" ledger E_B - E_A - Eg: {(d['EB']-d['EA']-d['Eg'])*MHz:+.3f} / {(o['EB']-o['EA']-o['Eg'])*MHz:+.3f} MHz")
