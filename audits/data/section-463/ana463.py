# Audit 463: state before the first photon vs the basin prediction, per run; populations vs P(l).
import sys, re, math
def parse(path):
    runs={}; cur=None
    for line in open(path):
        if line.startswith('START'):
            p=line.split(); cur=int(p[2]); runs[cur]=dict(n0=float(p[4]),L0=None,n=None,L=None,photon=False,end=False,closed=False)
        elif cur is None: continue
        elif line.startswith('WAVE_LOCK ') and not runs[cur]['photon']:
            nb=float(re.search(r'n_before=(\S+)',line).group(1)); na=float(re.search(r'n_after=(\S+)',line).group(1))
            Lh=float(re.search(r'L_hbar=(\S+)',line).group(1))
            r=runs[cur]
            if r['L0'] is None: r['L0']=Lh
            if r['n'] is not None and abs(nb-r['n'])>0.3: r['photon']=True; continue
            r['n']=na
        elif line.startswith('WAVE_LOCK_L') and not runs[cur]['photon']:
            runs[cur]['L']=float(re.search(r'L_after=(\S+)',line).group(1))
        elif line.startswith('END'):
            runs[cur]['end']=True; p=line.split(); runs[cur]['closed']=(p[4]=='4' and p[6]=='0')
    return runs
def predict(n0,L0):
    n=round(n0); return n, min(math.floor(L0), n-1)
out={}
for lab in sys.argv[1:]:
    runs=parse(f'audits/data/section-463/out/{lab}.err'); ok=0; tot=0; pop={}; res={}; nb=[0]
    for i,r in sorted(runs.items()):
        if not r['end'] or r['L'] is None: continue
        tot+=1; n,l=predict(r['n0'],r['L0']); nf=round(r['n']); lf=r['L']-0.5
        on=abs(r['n']-nf)<1e-6 and abs(lf-round(lf))<1e-6
        # engine closure before any photon (stop 4, 0 photons): the E1 channel is shut at the L rung -> s state
        if r['closed'] and r['L']<1.0: on=True; nf=n; lf=0.0
        good=on and nf==n and round(lf)==l; ok+=good
        basin=(nf==n and round(lf)==l)   # nearest rung, lock possibly unfinished
        nb[0]+=basin
        res[i]=(nf,round(lf))
        k=round(lf); pop[k]=pop.get(k,0)+1
        if not good: print(('  [unfinished lock]' if basin else '  [WRONG BASIN]'),f'  {lab} idx {i}: n0 {r["n0"]:.4f} L0 {r["L0"]:.4f} -> n {r["n"]:.6f} L {r["L"]:.6f} (predicted n {n} l {l})')
    nc=round(sum(r['n0'] for r in runs.values())/len(runs)); d=nc*nc+1/12
    P=[(2*l+1)/d for l in range(nc-1)]; P.append(1-sum(P))
    s='  '.join(f'l={l}: {pop.get(l,0)}/{tot} = {pop.get(l,0)/tot:.3f} (P {P[l]:.4f}, {(pop.get(l,0)/tot-P[l])/math.sqrt(P[l]*(1-P[l])/tot):+.1f} sigma; QM {(2*l+1)/nc**2:.3f})' for l in range(nc))
    print(f'{lab}: runs {tot}, on the predicted rung {ok}/{tot}, nearest rung = prediction {nb[0]}/{tot};  {s}')
    out[lab]=res
pairs=[k for k in out if k.startswith('H2-')]
for k in pairs:
    same=sum(1 for i in out[k] if i in out.get('H2',{}) and out[k][i]==out['H2'][i]); n=sum(1 for i in out[k] if i in out.get('H2',{}))
    print(f'paired {k} vs H2: same (n, l) {same}/{n}')
