# Audit 678: Ps 1S states at 1.0329 T from the engine's own prints (runs/ps) -- HFS from the field-independent sum and the Zeeman line f_Z,
# quantum (engine eigenstates) vs classical image (Ua + Za), channel by channel.   usage: python3 hfs678.py
import glob,collections,math
e=1.602176634e-19; h=6.62607015e-34; c=299792458.0; muB=9.2740100657e-24; alpha=7.2973525643e-3; g=2.00231930436092; B=1.0329
def three(f): return f'{f/1e9:.6f} GHz  {c/f*1e3:.6f} mm  {h*f/e*1e6:.6f} ueV'
S=collections.defaultdict(list); n=0
for f in sorted(glob.glob('runs/ps/*.err')):
    q=None
    for l in open(f):
        t=l.split()
        if l.startswith('BORN Q photon'):
            q=dict(E1S=float(t[t.index('E1S')+1]),Ua=float(t[t.index('Ua')+1]),Za=float(t[t.index('Za')+1]),Qa=float(t[t.index('Qa')+1]),m=int(t[t.index('m')+1]))
        elif l.startswith('CREM_STATE_ENERGY photon') and q is not None:
            q['Ra']=float(t[t.index('Ra')+1]); q['Rdef']=float(t[t.index('Rdef')+1])
        elif l.startswith('BORN field eigen (1S)') and q is not None:
            PS=float(l.split('P_S=')[1].split()[0]); n+=1
            S[(q['m'],round(q['E1S']/e*1e6,5),PS)].append(q); q=None
out=[f'{n} photons, {len(glob.glob("runs/ps/*.err"))} probes']
rows={}
for (m,E,PS),v in sorted(S.items()):
    sv=[-x['Qa']/e*1e6 for x in v]; rd=[x.get('Rdef',float('nan'))/e*1e6 for x in v]
    name='m=%+d'%m if m else ('m=0 p-type' if PS>=0.5 else 'm=0 o-type')
    rows[name]=(E,sum(sv)/len(sv))
    out.append(f'  {name:11s} P_S {PS:.6f} n {len(v):3d}  E1S {E:+.5f} ueV  S1S*=-Qa {sum(sv)/len(sv):+.5f} ueV (spread {max(sv)-min(sv):.1e})  Rdef {sum(rd)/len(rd):+.5f}  E_cl {E-sum(sv)/len(sv):+.5f}')
need=['m=+1','m=-1','m=0 o-type','m=0 p-type']
missing=[k for k in need if k not in rows]
if missing: out.append(f'  MISSING {missing}')
else:
    for lab,k in (('QM',lambda r:r[0]),('classical',lambda r:r[0]-r[1])):
        D=(k(rows['m=+1'])+k(rows['m=-1'])-k(rows['m=0 o-type'])-k(rows['m=0 p-type']))*1e-6*e/h
        fZ=(0.5*(k(rows['m=+1'])+k(rows['m=-1']))-k(rows['m=0 o-type']))*1e-6*e/h
        out.append(f'  {lab:9s} (E+1 + E-1) - (E0o + E0p) = {D/1e9:.6f} GHz ({h*D/e*1e6:.6f} ueV); - Ishida 203.3942: {(D/1e9-203.3942)*1e3:+.3f} MHz; - Ritter 203.38910: {(D/1e9-203.38910)*1e3:+.3f} MHz')
        out.append(f'            f_Z = E0o - E(+-1) = {three(-fZ)};  E+1 - E-1 = {(k(rows["m=+1"])-k(rows["m=-1"]))*1e-6*e/h/1e6:+.4f} MHz')
D=203.3942e9; qq=g*(1-5*alpha**2/24)*muB*B/(h*D); fZ=0.5*D*(math.sqrt(1+4*qq*qq)-1)
out.append(f'  Breit-Rabi with Ishida Delta_HFS and g\' at {B} T (measurement formula, not a measurement at this field): f_Z = {three(fZ)}')
print('\n'.join(out))
