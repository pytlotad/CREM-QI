# Audit 475: w at preparation, before and after the 2p -> 1s photon, and the lifetime spectrum it implies.
import math, glob
def unit(v): n=math.sqrt(sum(x*x for x in v)); return [x/n for x in v]
def wof(m1,m2):  # w = (1 - s1.s2)/2, s1 = -m1^ (electron), s2 = +m2^ (positron)
    a=unit(m1); b=unit(m2); return (1+sum(x*y for x,y in zip(a,b)))/2
starts=set(); photons=[]
for f in sorted(glob.glob('run/err*.txt')):
    for line in open(f):
        t=line.split()
        if line.startswith('CREM_LEDGER_START'): starts.add(tuple(t[2:5]+t[6:9]))
        elif line.startswith('CREM_LEDGER '):
            bars=[k for k,x in enumerate(t) if x=='|']
            def side(s):
                seg=t[s:]; i1=seg.index('m1'); i2=seg.index('m2'); return wof([float(x) for x in seg[i1+1:i1+4]],[float(x) for x in seg[i2+1:i2+4]])
            photons.append((side(bars[0]),side(bars[1]),float(t[t.index('Eg')+1])))
w0=[wof([float(x) for x in st[:3]],[float(x) for x in st[3:]]) for st in starts]
print(f'unique start states {len(w0)}, photons {len(photons)}')
dw=max((abs(p[1]-p[0]) for p in photons),default=0)
print(f'H475a: max |w_after - w_before| at the photon = {dw:.3e}')
wb=[p[0] for p in photons]; n=len(wb)
hist=lambda ws:[sum(1 for w in ws if i/10<=w<(i+1)/10 or (i==9 and w==1)) for i in range(10)]
edge=sum(1 for w in wb if w<0.05 or w>0.95)/n; edge0=sum(1 for w in w0 if w<0.05 or w>0.95)/len(w0)
print(f'H475b: fraction with w < 0.05 or > 0.95: start {edge0:.3f} (n {len(w0)}), before the photon {edge:.3f} (n {n}); uniform 0.10, sigma {math.sqrt(0.09/n):.3f}')
print('       histogram w at start:',hist(w0)); print('       histogram w before photon:',hist(wb))
done=[None]*n
tp=124.494e-12; to=138.674e-9
wf=[p[1] for p in photons]
taus=[1/(w/tp+(1-w)/to) for w in wf]
long=sum(1 for x in taus if x>10e-9)/n
print(f'H475c: fraction of pairs with lifetime > 10 ns: {long:.3f} (QM, uncorrelated formation: 0.75; measured intensities 1 : 3)')
print('       lifetimes [ns] sorted (quantiles 10/25/50/75/90 %):',[round(sorted(taus)[int(q*(n-1))]*1e9,4) for q in (0.1,0.25,0.5,0.75,0.9)])
for t in (0.5e-9,2e-9,10e-9,50e-9,142e-9):
    S=sum(math.exp(-t/x) for x in taus)/n; Sqm=sum((w/2)*math.exp(-t/tp)+(1-w/2)*math.exp(-t/to) for w in wf)/n
    print(f'       survival S({t*1e9:6.1f} ns): model {S:.4f}   QM (P_S = w/2, two components) {Sqm:.4f}')
