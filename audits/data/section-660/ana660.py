# Audit 660: antihydrogen 2P -> 1S lines at 1.0329 T against ALPHA (Nature 578, 375 (2020)).
# usage: python3 ana660.py <runs> [<runs without the new default>]
import re,sys,math,os
h=6.62607015e-34; hbar=h/(2*math.pi)
F=r'(-?[0-9.e+-]+)'
def parse(f):
    if not os.path.exists(f): return None
    t=open(f).read()
    b=re.search(r'BHAT '+F+' '+F+' '+F,t); m=re.search(r'CREM_LEDGER t \S+ Eg (\S+) .*? sphoton '+F+' '+F+' '+F+r' .*\| A .*? m1 '+F+' '+F+' '+F+' m2 '+F+' '+F+' '+F,t)
    if not (b and m): return None
    B=[float(x) for x in b.groups()]; v=[float(x) for x in m.groups()]
    Eg=v[0]; sp=v[1:4]; m1=v[4:7]; m2=v[7:10]
    dot=lambda a,c: sum(x*y for x,y in zip(a,c)); nrm=lambda a: math.sqrt(dot(a,a))
    dm=dot(sp,B)/hbar
    sA=-dot(m1,B)/nrm(m1)   # antiproton: spin antiparallel to its moment (q < 0, g > 0)
    sP=+dot(m2,B)/nrm(m2)   # positron: spin parallel to its moment
    return dict(Eg=Eg,dm=dm,sA=sA,sP=sP)
def report(d,tag):
    R=[parse(f'{d}/A{ph}_{i}.err') for ph in (1,2) for i in range(64)]
    R=[r for r in R if r]
    print(f'=== {tag}: photons {len(R)}/128')
    from collections import Counter
    print('  (Delta m, positron spin sign, antiproton spin sign) counts:',dict(Counter((round(r["dm"]),r["sP"]>0,r["sA"]>0) for r in R)))
    S=sorted(r['Eg']/h/1e6 for r in R if round(r['dm'])==0 and r['sP']<0 and r['sA']<0)
    print(f'  selected (pi, positron and antiproton spins antiparallel to B): {len(S)}')
    if len(S)<2: return
    print('  lines (MHz): '+' '.join(f'{x:.1f}' for x in S))
    gaps=[(S[i+1]-S[i],i) for i in range(len(S)-1)]; g,i=max(gaps)
    lo,hi=S[:i+1],S[i+1:]
    mean=lambda x: sum(x)/len(x); rms=lambda x: math.sqrt(sum((y-mean(x))**2 for y in x)/len(x))
    spread=max(rms(lo),rms(hi))
    sep=mean(hi)-mean(lo)
    print(f'  largest gap {g:.1f} MHz between groups of {len(lo)} and {len(hi)}; rms {rms(lo):.1f} / {rms(hi):.1f} MHz; gap/spread {g/max(spread,1e-9):.1f}')
    print(f'  H660b: separation of the group centres {sep/1e3:.4f} GHz (measured 14.945(75) GHz; window 14.795 .. 15.095); two groups (gap > 3 spread): {g>3*spread}')
    print(f'  H660c: centres {mean(hi):.1f} / {mean(lo):.1f} MHz against 2466051189(76) / 2466036395(81): {mean(hi)-2466051189:+.1f} / {mean(lo)-2466036395:+.1f} MHz')
report(sys.argv[1],'default (CREM_PAIR_TRANSPORT)')
if len(sys.argv)>2: report(sys.argv[2],'CREM_NO_PAIR_TRANSPORT=1 (H660d)')
