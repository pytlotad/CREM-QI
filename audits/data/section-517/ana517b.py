# Audit 517: production photons, default (k = 1 in the S band) vs CREM_NO_SOMMERFELD_S=1 (import 17); H and Ps, start n = 2 with LFRAC 0.75 (2p, L' = hbar/2) and 0.65 (L = 1.3 hbar, L' = 0.3 hbar), and
# LFRAC 0.25 (L = hbar/2, ladder) and 0.4 (L = 0.8 hbar, off the ladder).  Predicted Eg(new) - Eg(old) = dE_2(L) - dE_1(L'),
# dE_n(L) = B_n (1 - 1/(2L)), B_n from ana517.txt (bare E_rel(1/2) - E_rel(1)).  python3 ana517b.py
import math
h=6.62607015e-34; MHz=1e-6/h; hb=1.054571817e-34
B={'pr':{2:-21886.648,1:-175093.183},'el':{2:-10949.284,1:-87594.271}}
def rows(f):
    out=[]
    for line in open(f):
        if not line.startswith('CREM_LEDGER '): continue
        t=line.split(); bars=[k for k,x in enumerate(t) if x=='|']; head=t[:bars[0]]; Bf=t[bars[0]+1:bars[1]]; A=t[bars[1]+1:]
        g=lambda seg,key: float(seg[seg.index(key)+1])
        out.append(dict(t=g(head,'t'),Eg=g(head,'Eg'),LB=g(Bf,'L')/hb,LA=g(A,'L')/hb,UFB=g(Bf,'UF'),UFA=g(A,'UF')))
    return out
dE=lambda p,n,L: B[p][n]*(1-1/(2*L)) if L<1 else 0.0
for p,lab in (('pr','H'),('el','Ps')):
  for lf in ('0.25','0.4','0.75','0.65'):
    N=rows(f'prod/new_{p}_{lf}.err'); O=rows(f'prod/old_{p}_{lf}.err')
    print(f"{lab} LFRAC {lf}: photons new {len(N)} old {len(O)}")
    for i,(a,b) in enumerate(zip(N,O)):
      pred=dE(p,2,a['LB'])-dE(p,1,a['LA'])
      print(f"  {i}: L {a['LB']:.6f} -> {a['LA']:.6f} hbar (old {b['LB']:.6f} -> {b['LA']:.6f}); t {a['t']*1e9:.5f} / {b['t']*1e9:.5f} ns;"
            f" UF_B {a['UFB']*MHz:+.2f} / {b['UFB']*MHz:+.2f}; Eg new - old {(a['Eg']-b['Eg'])*MHz:+.4f} MHz, predicted {pred:+.4f}, diff {(a['Eg']-b['Eg'])*MHz-pred:+.2e}")
