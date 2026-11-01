# Audit 677: 1S positron spin-flip lines of antihydrogen at 1.0329 T from the engine's own prints (../section-676/runs) vs ALPHA 2017.
# usage: python3 hfs677.py   (run from this directory)
import glob,re,collections
e=1.602176634e-19; h=6.62607015e-34; c=299792458.0
def states(d):
    S=collections.defaultdict(list)
    for f in sorted(glob.glob(f'../section-676/runs/{d}/*.err')):
        E1S=PS=None
        for l in open(f):
            if l.startswith('PAIR eigen draw'):
                t=l.split(); E1S=float(t[t.index('E1S')+1]); PS=float(t[t.index('P_S')+1])
            if l.startswith('PAIR eigen ledger') and E1S is not None:
                t=l.split(); s=float(t[t.index('S1S')+1])-float(t[t.index('Zsa')+1])
                S[(round(E1S/e*1e6,6),round(PS,6))].append(s)
    return S
def show(f): return f'{f/1e9:.6f} GHz  {c/f*1e3:.6f} mm  {h*f/e*1e6:.6f} ueV'
out=[]
for d in ('on','off'):
    S=states(d); out.append(f'== {d}/ : {sum(len(v) for v in S.values())} photons')
    rows=[]
    for (E,PS),v in sorted(S.items()):
        sp=max(v)-min(v); s=sum(v)/len(v)/e*1e6
        rows.append((E,PS,s)); out.append(f'  E1S {E:+.6f} ueV  P_S {PS:.6f}  n {len(v):3d}  S1S* {s:+.6f} ueV (spread {sp/e*1e6:.1e})  E_cl {E-s:+.6f} ueV')
    pure=[r for r in rows if r[1]==0.0]; mix=[r for r in rows if r[1]>0]
    assert len(pure)==2 and len(mix)==2, rows
    for lab,k in (('QM',lambda r:r[0]),('classical',lambda r:r[0]-r[2])):
        D=(sum(k(r) for r in pure)-sum(k(r) for r in mix))*1e-6*e/h
        out.append(f'  {lab:9s} (E_b+E_d)-(E_a+E_c) = {D/1e6:.4f} MHz  ({h*D/e*1e6:.6f} ueV)  - ALPHA 1420.4(5): {D/1e6-1420.4:+.4f} MHz = {(D/1e6-1420.4)/0.5:+.1f} sigma; - H 1420.405752: {D/1e6-1420.405752:+.4f} MHz')
        # ALPHA labels: in antihydrogen the positron moment is along its spin, so |up> lies low; a/b = positron up (lower pair), c/d = positron down.
        a,b,cc,dd=sorted(pure+mix,key=lambda r:r[0])  # labels from the QM order a<b<c<d at high field (ALPHA Fig. 2), kept for the classical image
        fcb=(k(cc)-k(b))*1e-6*e/h; fda=(k(dd)-k(a))*1e-6*e/h
        out.append(f'    f_cb {show(fcb)}   (c P_S {cc[1]:.4f}, b P_S {b[1]:.4f})')
        out.append(f'    f_da {show(fda)}   (d P_S {dd[1]:.4f}, a P_S {a[1]:.4f})')
print('\n'.join(out))
