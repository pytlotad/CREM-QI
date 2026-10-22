# Audit 557: model cascades (CREM_LEDGER lines) -> per-capture line probabilities, compared with cascade557 (QM).  python3 ana557.py dir samples-file label
import sys, os, math, re
d,samples,label=sys.argv[1],sys.argv[2],sys.argv[3]
B=13.5983*1.602176634e-19
S=[l.split() for l in open(samples)]
N=len(S); cnt={'direct1s':0,'lya':0,'end2s':0,'lyb+':0,'ha':0,'hb':0,'pa':0,'endother':0,'unfinished':0}; times=[]; nph=[]; stuck=[]
for i,n,L in S:
    n=int(n)
    if n==1: cnt['direct1s']+=1; continue
    f=os.path.join(d,f'c{i}.err')
    if not os.path.exists(f): cnt['unfinished']+=1; continue
    txt=open(f).read()
    if 'END idx' not in txt or 'stop 0 ' in txt: cnt['unfinished']+=1; continue   # stop 0 = censored by the budget, not a result
    last=n; ph=0; t=0.0
    for l in txt.splitlines():
        if not l.startswith('CREM_LEDGER t'): continue
        tok=l.split(); bars=[k for k,x in enumerate(tok) if x=='|']
        Eb=float(tok[tok.index('E',bars[0])+1]); Ea=float(tok[tok.index('E',bars[1])+1]); t=float(tok[2])
        nb=round(math.sqrt(B/abs(Eb))); na=round(math.sqrt(B/abs(Ea))); ph+=1
        if (nb,na)==(2,1): cnt['lya']+=1
        elif na==1: cnt['lyb+']+=1
        if (nb,na)==(3,2): cnt['ha']+=1
        if (nb,na)==(4,2): cnt['hb']+=1
        if (nb,na)==(4,3): cnt['pa']+=1
        last=na
    if last==2: cnt['end2s']+=1
    elif last!=1:
        cnt['endother']+=1; stuck.append((i,n,L,last,re.search(r'stop (\d+)',txt).group(1)))
    times.append(t); nph.append(ph)
done=N-cnt['unfinished']
print(f"{label}: samples {N}, finished {done}")
print(f"  per capture: direct 1s {cnt['direct1s']/done:.3f}; Lyman-alpha {cnt['lya']/done:.3f}; ends in 2s {cnt['end2s']/done:.3f}; Lyman-beta+ {cnt['lyb+']/done:.3f}; stops above n = 2 {cnt['endother']/done:.3f}")
print(f"  Halpha {cnt['ha']/done:.3f}, Hbeta {cnt['hb']/done:.3f}, Halpha/Hbeta {cnt['ha']/max(cnt['hb'],1):.2f}; Paschen-alpha {cnt['pa']/done:.3f}; mean photons (n >= 2) {sum(nph)/max(len(nph),1):.2f}; median last-photon time {sorted(times)[len(times)//2]*1e9 if times else float('nan'):.2f} ns")
print("  stopped above n = 2 (index, start n, start L, last n, stop cause): "+str(stuck))
