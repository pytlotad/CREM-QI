# Audit 355: singlet transport -- w, L and QR lifetime at n = 1 after a Langer 2p start; compare o-Ps END with 353.
import glob,os,sys,statistics as st
D=sys.argv[1] if len(sys.argv)>1 else '.'; R353=sys.argv[2] if len(sys.argv)>2 else None
T0=124.4942e-12; EPS=1/1113.8978
out=open(os.path.join(D,'analysis.txt'),'w')
def P(*a): print(*a); print(*a,file=out)
rows=[]
for f in sorted(glob.glob(os.path.join(D,'out','*.txt'))):
    for l in open(f):
        if l.startswith('END'):
            x=l.split(); d=dict(t.split('=') for t in x[4:]); c=float(d['n/psi0']); G=float(d['Gf']); br=G*T0/c
            n=1/float(d['nE']); L=float(d['L'])
            rows.append(dict(i=int(x[2]),ch='p-Ps' if x[3]=='para' else 'o-Ps',stop=d['stop'],n=n,L=L,w=(br-EPS)/(1-EPS),tau=T0*2*n**3*L/br,line=l.strip()))
for ch in ['p-Ps','o-Ps']:
    R=[r for r in rows if r['ch']==ch]
    t=[r['tau'] for r in R]
    P(f"{ch}: N={len(R)} stops {[r['stop'] for r in R]}; n {sorted({round(r['n'],6) for r in R})}; max|L-0.5| {max(abs(r['L']-0.5) for r in R):.2e}")
    P(f"   w min {min(r['w'] for r in R):.6f} max {max(r['w'] for r in R):.6f}; tau mean {st.mean(t):.5e} [{min(t):.5e}..{max(t):.5e}]; <1/G><G> {st.mean(t)*st.mean(1/x for x in t):.6f}")
tp=st.mean(r['tau'] for r in rows if r['ch']=='p-Ps'); to=st.mean(r['tau'] for r in rows if r['ch']=='o-Ps')
P(f"tau_p/125.14ps = {tp/125.14e-12:.4f}; tau_o/142.04ns = {to/142.04e-9:.4f}; tau_o/tau_p = {to/tp:.2f} (1/eps 1113.90, measured 1135.0)")
if R353:
    old={l.split()[2]:l.strip() for f in glob.glob(os.path.join(R353,'out','n2-j0.75-A-*.txt')) for l in open(f) if l.startswith('END') and 'ortho' in l}
    new={r['line'].split()[2]:r['line'] for r in rows if r['ch']=='o-Ps'}
    P(f"o-Ps END identical to 353 (A, 2p): {sum(old.get(k)==v for k,v in new.items())}/{len(new)}")
h=os.path.join(D,'hold0.txt')
if os.path.exists(h):
    ws={}; cur=None
    for l in open(h):
        if l.startswith('BEGIN'): cur=l.split()[2]; ws[cur]=[]
        elif l.startswith('CREM_HOLD'): x=l.split(); ws[cur].append((float(x[1]),float(x[2]),float(x[4])))
    for k,v in ws.items():
        if v: P(f"hold idx 0 {'p-Ps' if k=='1' else 'o-Ps'}: {len(v)} samples over {v[-1][0]*1e12:.1f} ps; w min {min(x[2] for x in v):.6f} max {max(x[2] for x in v):.6f}; L {min(x[1] for x in v):.5f}..{max(x[1] for x in v):.5f}")
