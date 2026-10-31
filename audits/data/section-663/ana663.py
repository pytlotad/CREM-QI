# Audit 663: natural width from the engine's own emission rate (CREM_NATURAL_WIDTH=1).
# usage: python3 ana663.py <dir with runs/ b0/ ps/ eq/> <662 dir>
import re,sys,math,os
h=6.62607015e-34; c0=299792458.0; e0=1.602176634e-19
line=lambda f: f'{f:.1f} MHz ({c0/(f*1e6)*1e9:.7f} nm; {h*f*1e6/e0:.9f} eV)'
ueV=lambda f: f'{h*f*1e6/e0*1e6:.4f} ueV'
d,old=sys.argv[1],sys.argv[2]
A_H,A_Ps=6.2649e8,3.1342e8          # exact QM E1 2P -> 1S with reduced mass (prereg)
def rd(f):
    if not os.path.exists(f): return None
    t=open(f).read()
    w=re.search(r'NATURAL width Gamma (\S+) hbarGamma (\S+) delta (\S+)',t); g=re.search(r'CREM_LEDGER t \S+ Eg (\S+)',t)
    lab=re.findall(r'^PAIR (?:eigen draw|field eigen).*$',t,re.M)
    return dict(G=float(w.group(1)) if w else None,hG=float(w.group(2)) if w else None,dl=float(w.group(3)) if w else 0.0,
                Eg=float(g.group(1)) if g else None,lab=lab)
def binom(k,n,p): n=max(n,1); s=math.sqrt(p*(1-p)/n); return k/n,s,abs(k/n-p)<=2*s
X=[]
for sub in ('runs','b0'):
    R=[r for r in (rd(f'{d}/{sub}/A1_{i}.err') for i in range(64)) if r and r['G']]
    X+=[r['dl']/r['hG'] for r in R]
    G=[r['G'] for r in R]
    if G:
        m=sum(G)/len(G); print(f'=== antihydrogen {sub}: {len(G)} photons, Gamma engine {m:.5e} 1/s (min {min(G):.5e}, max {max(G):.5e}), '
            f'Gamma/2pi {m/2/math.pi/1e6:.2f} MHz ({ueV(m/2/math.pi/1e6)}), tau {1e9/m:.4f} ns; QM {A_H:.4e}: {100*(m/A_H-1):+.2f} % -> H663b {abs(m/A_H-1)<0.01}')
    # H663d: labels and Eg - delta against 662, same seed
    same=0; eg=0; n=0
    for i in range(64):
        a=rd(f'{d}/{sub}/A1_{i}.err'); b=rd(f'{old}/{sub}/A1_{i}.err')
        if not (a and b and a['Eg'] and b['Eg']): continue
        n+=1; same+=a['lab']==b['lab']; eg+=abs(a['Eg']-a['dl']-b['Eg'])<=1e-12*b['Eg']
    print(f'  H663d {sub}: labels = 662 {same}/{n}, Eg - delta = 662 Eg {eg}/{n}')
    # line shapes: group photons by the 662 sharp line of the same seed
    grp={}
    for i in range(64):
        a=rd(f'{d}/{sub}/A1_{i}.err'); b=rd(f'{old}/{sub}/A1_{i}.err')
        if a and b and a['Eg'] and b['Eg']: grp.setdefault(round(b['Eg']/h/1e6,1),[]).append((a['Eg']-b['Eg'])/h/1e6)
    for v in sorted(grp):
        s=sorted(grp[v]); med=s[len(s)//2]
        q=(s[(3*len(s))//4]-s[len(s)//4]) if len(s)>=4 else float('nan')
        print(f'    {line(v)}: {len(s)} photons, median shift {med:+.1f} MHz ({ueV(med)}), IQR {q:.1f} MHz'+('' if q!=q else f' ({ueV(q)})'))
# post hoc (uzupelnienie 2): the stream splitMix64(seed ^ 0x663 ^ photon) gives the SAME u for the same seed in runs/ and b0/ (and ps/),
# so only one x per seed is independent: H663a over the seeds, x from b0/ (or runs/ where b0/ has none).
U={}
for sub in ('b0','runs'):
    for i in range(64):
        r=rd(f'{d}/{sub}/A1_{i}.err')
        if r and r['G'] and i not in U: U[i]=r['dl']/r['hG']
dup=sum(1 for i in range(64) for a,b in [(rd(f'{d}/runs/A1_{i}.err'),rd(f'{d}/b0/A1_{i}.err'))] if a and b and a['G'] and b['G'] and abs(a['dl']/a['hG']-b['dl']/b['hG'])<1e-9)
k=sum(1 for x in U.values() if abs(x)<0.5); f,s,ok=binom(k,len(U),0.5)
print(f'=== H663a: |delta|/(hbar Gamma) < 1/2 in {k}/{len(U)} independent seeds = {f:.3f} (Cauchy 0.5 +- {2*s:.3f}): {ok}'
      f'  [all {len(X)} photons: {sum(1 for x in X if abs(x)<0.5)}; same x in runs/ and b0/ for {dup} seeds]')
for ph,tag in ((1,'p-Ps'),(2,'o-Ps')):
    R=[r for r in (rd(f'{d}/ps/P{ph}_{i}.err') for i in range(16)) if r and r['G']]
    if not R: print(f'=== {tag}: no photons'); continue
    G=[r['G'] for r in R]; m=sum(G)/len(G)
    print(f'=== {tag} (PHEN {ph}, 1e-6 T): {len(G)} photons, Gamma engine {m:.5e} 1/s (min {min(G):.5e}, max {max(G):.5e}), Gamma/2pi {m/2/math.pi/1e6:.2f} MHz '
          f'({ueV(m/2/math.pi/1e6)}), tau {1e9/m:.4f} ns; QM {A_Ps:.4e}: {100*(m/A_Ps-1):+.2f} % -> H663c {abs(m/A_Ps-1)<0.01}')
eq=sum(1 for i in range(4) if os.path.exists(f'{d}/eq/A1_{i}.err') and open(f'{d}/eq/A1_{i}.err','rb').read()==open(f'{old}/b0/A1_{i}.err','rb').read())
print(f'=== eq: without CREM_NATURAL_WIDTH = 662 b0 byte for byte {eq}/4')
