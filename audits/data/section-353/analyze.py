# Audit 353: Langer starts L = (l + 1/2) hbar; per config, variant, channel: photon path, terminal n, L, QR lifetime.
import glob,os,sys,math,statistics as st
D=sys.argv[1] if len(sys.argv)>1 else '.'
T0=124.4942e-12; EPS=1/1113.8978; eV=1.602176634e-19
rows=[]
for f in sorted(glob.glob(os.path.join(D,'out','*.txt'))):
    nl,jj,v,i=os.path.basename(f)[:-4].split('-'); n0=int(nl[1:]); j0=float(jj[1:]); L0=j0*n0
    cur=None
    for l in open(f):
        if l.startswith('BEGIN'): cur=dict(n0=n0,L0=round(L0,3),v=v,i=int(i),ch='p-Ps' if l.split()[2]=='1' else 'o-Ps',ph=[])
        elif l.startswith('CREM_REACH'):
            x=l.split()[1:]; cur['ph'].append((float(x[4])/eV,int(x[10]),float(x[11]),float(x[12])*float(x[1])))
        elif l.startswith('END'):
            d=dict(t.split('=') for t in l.split()[4:]); cur['stop']=d['stop']; cur['n']=1/float(d['nE']); cur['L']=float(d['L'])
            c=float(d['n/psi0']); G=float(d['Gf']); br=G*T0/c if c>0 else float('nan')
            cur['w']=(br-EPS)/(1-EPS); cur['tau']=T0*2*cur['n']**3*cur['L']/br if br>0 else float('nan'); rows.append(cur)
out=open(os.path.join(D,'analysis.txt'),'w')
def P(*a): print(*a); print(*a,file=out)
P("n0 L0   var ch   N  stops                 photons(E eV,k)            terminal n         terminal L (min..max)    tau QR mean [min..max]   <1/G><G>")
keys=sorted({(r['n0'],r['L0'],r['v'],r['ch']) for r in rows})
for k in keys:
    R=[r for r in rows if (r['n0'],r['L0'],r['v'],r['ch'])==k]
    stops={}; [stops.__setitem__(r['stop'],stops.get(r['stop'],0)+1) for r in R]
    paths={}; [paths.__setitem__(' '.join(f"{e:.3f}({kk})" for e,kk,_,_ in r['ph']) or '-',paths.get(' '.join(f"{e:.3f}({kk})" for e,kk,_,_ in r['ph']) or '-',0)+1) for r in R]
    ns=sorted({round(r['n'],5) for r in R}); Ls=[r['L'] for r in R]
    t=[r['tau'] for r in R if r['stop']!='barrier' and math.isfinite(r['tau'])]
    tm=st.mean(t) if t else float('nan'); spread=(st.mean(t)*st.mean(1/x for x in t)) if t else float('nan')
    P(f"{k[0]}  {k[1]:.3f} {k[2]}  {k[3]} {len(R)}  {stops}  {paths}  n={ns}  L {min(Ls):.4f}..{max(Ls):.4f}  tau {tm:.4e} [{min(t) if t else float('nan'):.3e}..{max(t) if t else float('nan'):.3e}]  {spread:.5f}")
P("\nAt n = 1 (not barrier): |L - 1/2| max, tau by channel, w")
for ch in ['p-Ps','o-Ps']:
    R=[r for r in rows if r['ch']==ch and abs(r['n']-1)<1e-6 and r['stop']!='barrier']
    if not R: continue
    t=[r['tau'] for r in R]
    P(f"{ch}: N={len(R)} max|L-0.5| {max(abs(r['L']-0.5) for r in R):.2e}; tau mean {st.mean(t):.4e} s median {st.median(t):.4e} [{min(t):.4e}..{max(t):.4e}]; <1/G><G> {st.mean(t)*st.mean(1/x for x in t):.5f}; w min {min(r['w'] for r in R):.4g} mean {st.mean(r['w'] for r in R):.4g}")
