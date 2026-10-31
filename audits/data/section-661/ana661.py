# Audit 661: antihydrogen 2P -> 1S sigma lines from 1Sd at 1.0329 T with CREM_PAIR_FIELD_EIGEN=1 against ALPHA (Nature 578, 375 (2020)).
# usage: python3 ana661.py <runs> [<660 runs (default engine, no flag)>]
import re,sys,math,os
h=6.62607015e-34; hbar=h/(2*math.pi); F=r'(-?[0-9.e+-]+)'
ZERO=8176.0   # 660d: model minus measurement without field (1S Lamb shift), MHz
ALPHA=(2466051189.0,76.0,2466036395.0,81.0)
c0=299792458.0; e0=1.602176634e-19
line=lambda f: f'{f:.1f} MHz ({c0/(f*1e6)*1e9:.7f} nm; {h*f*1e6/e0:.9f} eV)'   # absolute line: vacuum wavelength, photon energy
shift=lambda d: f'{d:+.1f} MHz ({h*d*1e6/e0*1e6:+.4f} ueV)'                     # difference
def parse(f):
    if not os.path.exists(f): return None
    t=open(f).read()
    b=re.search(r'BHAT '+F+' '+F+' '+F,t); m=re.search(r'CREM_LEDGER t \S+ Eg (\S+) .*? sphoton '+F+' '+F+' '+F+r' .*\| A .*? m1 '+F+' '+F+' '+F+' m2 '+F+' '+F+' '+F,t)
    if not (b and m): return None
    B=[float(x) for x in b.groups()]; v=[float(x) for x in m.groups()]
    dot=lambda a,c: sum(x*y for x,y in zip(a,c)); nrm=lambda a: math.sqrt(dot(a,a))
    r=dict(Eg=v[0]/h/1e6,dm=round(dot(v[1:4],B)/hbar),sA=-dot(v[4:7],B)/nrm(v[4:7]),sP=dot(v[7:10],B)/nrm(v[7:10]))
    p=re.search(r'PAIR eigen photon k=(\d+) .* import '+F,t); r['k'],r['imp']=(int(p.group(1)),float(p.group(2))/h/1e6) if p else (None,None)
    return r
mean=lambda x: sum(x)/len(x); rms=lambda x: math.sqrt(sum((y-mean(x))**2 for y in x)/len(x)) if x else float('nan')
def report(d,tag):
    R=[r for r in (parse(f'{d}/A{ph}_{i}.err') for ph in (1,2) for i in range(64)) if r]
    print(f'=== {tag}: photons {len(R)}/128, eigen path {sum(1 for r in R if r["k"] is not None)}')
    imp=[r['imp'] for r in R if r['imp'] is not None]
    if imp: print(f'  import energy (eigen - classical) MHz: rms {math.sqrt(mean([x*x for x in imp])):.1f}, max |.| {max(abs(x) for x in imp):.1f}')
    S=sorted(r['Eg'] for r in R if r['dm']==1 and r['sP']<0 and r['sA']<0)
    print(f'  sigma Delta m=+1 from 1Sd: {len(S)} lines: '+'; '.join(line(x) for x in sorted(set(round(y,1) for y in S))))
    if len(S)<2: return
    g,i=max((S[k+1]-S[k],k) for k in range(len(S)-1)); lo,hi=S[:i+1],S[i+1:]; sp=max(rms(lo),rms(hi))
    two=g>3*sp and sp<100.0; sep=(mean(hi)-mean(lo))/1e3
    print(f'  largest gap {g:.1f} MHz, groups {len(lo)} / {len(hi)}, rms {rms(lo):.3f} / {rms(hi):.3f} MHz')
    print(f'  H661b: two groups {two}; separation {sep:.4f} GHz = {h*sep*1e9/e0*1e6:.3f} ueV (14.945(75), window 14.795 .. 15.095): {two and 14.795<=sep<=15.095}')
    c,f=mean(hi)-ZERO,mean(lo)-ZERO
    print(f'  H661c: centres minus model zero {line(c)} / {line(f)}; minus ALPHA {shift(c-ALPHA[0])} / {shift(f-ALPHA[2])} (window 300 MHz): {abs(c-ALPHA[0])<=300 and abs(f-ALPHA[2])<=300}')
    print(f'  H661d: 2Pf- share {len(lo)}/{len(S)} = {len(lo)/len(S):.3f} (uniform-projection theory 0.14)')
    from collections import Counter
    print('  all lines by (dm, e+ spin up, pbar spin up):')
    for key,n in sorted(Counter((r['dm'],r['sP']>0,r['sA']>0) for r in R).items()):
        L=sorted(r['Eg'] for r in R if (r['dm'],r['sP']>0,r['sA']>0)==key)
        vals=sorted(set(round(x,1) for x in L))
        print(f'    {key}: {n} photons, distinct lines {len(vals)}: '+' '.join(f'{x:.1f}' for x in vals[:8]))
report(sys.argv[1],'CREM_PAIR_FIELD_EIGEN=1')
if len(sys.argv)>2: report(sys.argv[2],'660 default (no flag)')
