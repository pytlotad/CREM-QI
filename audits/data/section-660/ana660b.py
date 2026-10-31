# Audit 660, POST HOC (uzupelnienie 1): the ALPHA lines 1Sd -> 2Pc- / 2Pf- are sigma lines, not pi (theory660.py); select Delta m = +1
# (photon takes +hbar along B: positron m_l +1 -> 0, spin conserved), positron and antiproton spins antiparallel to B.  Also report the
# pi lines from 1Sd (theory separation 20.85 GHz) and the zero-field lines.  usage: python3 ana660b.py <runs> <runs0> <zero-field runs>
import sys,math
sys.argv=[sys.argv[0]]+sys.argv[1:]
from importlib.machinery import SourceFileLoader
import os
h=6.62607015e-34
src=open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'ana660.py')).read().split('def report')[0]
exec(src)
mean=lambda x: sum(x)/len(x); rms=lambda x: math.sqrt(sum((y-mean(x))**2 for y in x)/len(x))
def lines(d,dm):
    R=[parse(f'{d}/A{ph}_{i}.err') for ph in (1,2) for i in range(64)]
    return sorted(r['Eg']/h/1e6 for r in R if r and round(r['dm'])==dm and r['sP']<0 and r['sA']<0)
def split(S,ref,tag):
    print(f'  {tag}: {len(S)} lines (MHz): '+' '.join(f'{x:.1f}' for x in S))
    if len(S)<2: return
    g,i=max((S[k+1]-S[k],k) for k in range(len(S)-1)); lo,hi=S[:i+1],S[i+1:]; sp=max(rms(lo),rms(hi))
    print(f'    largest gap {g:.1f} MHz, groups {len(lo)} / {len(hi)}, rms {rms(lo):.1f} / {rms(hi):.1f} MHz, two groups {g>3*sp}; separation {(mean(hi)-mean(lo))/1e3:.4f} GHz (ref {ref})')
    print(f'    centres {mean(hi):.1f} / {mean(lo):.1f} MHz; minus ALPHA 2466051189 / 2466036395: {mean(hi)-2466051189:+.1f} / {mean(lo)-2466036395:+.1f} MHz')
for d,tag in ((sys.argv[1],'default'),(sys.argv[2],'CREM_NO_PAIR_TRANSPORT=1')):
    print(f'=== {tag}')
    split(lines(d,+1),'measured 14.945(75) GHz, window 14.795 .. 15.095','sigma Delta m=+1 (ALPHA 2Pc-/2Pf-)')
    split(lines(d,0),'theory 20.85 GHz (theory660.py)','pi Delta m=0 from 1Sd')
if len(sys.argv)>3:
    print('=== zero field (BEXT=1e-6, 16 probes, default)')
    Z=[]
    for ph in (1,2):
        for i in range(8):
            t=open(f'{sys.argv[3]}/A{ph}_{i}.err').read(); m=re.search(r'CREM_LEDGER t \S+ Eg (\S+)',t)
            if m: Z.append((float(m.group(1))/h/1e6,f'A{ph}_{i}'))
    for v,n in sorted(Z): print(f'  {n:6s} {v:.1f} MHz')
