# Audit 351: |L| and w at n = 1 during a 0.5 ns hold (no emission), per trajectory and channel.
import glob,os,sys,statistics as st
D=sys.argv[1] if len(sys.argv)>1 else '.'
out=open(os.path.join(D,'analysis.txt'),'w')
def P(*a): print(*a); print(*a,file=out)
res=[]
for f in sorted(glob.glob(os.path.join(D,'out','*.txt')),key=lambda x:int(os.path.basename(x)[:-4])):
    cur=None
    for l in open(f):
        if l.startswith('BEGIN'):
            _,i,ph=l.split(); cur=dict(i=int(i),ch='p-Ps' if ph=='1' else 'o-Ps',t=[],L=[],n=[],w=[]); res.append(cur)
        elif l.startswith('CREM_HOLD'):
            x=list(map(float,l.split()[1:])); cur['t'].append(x[0]); cur['L'].append(x[1]); cur['n'].append(x[2]); cur['w'].append(x[3])
        elif l.startswith('END'): cur['end']=l.split()[4]
P("idx ch   N    span(ps)  L0      max|dL|   slope(hbar/ns)  t(0.5hbar)/ns  n range           w min    w mean   w max")
for r in res:
    if len(r['t'])<3: P(f"{r['i']:3d} {r['ch']} too few samples ({len(r['t'])})"); continue
    t=r['t']; L=r['L']; tm=st.mean(t); Lm=st.mean(L)
    slope=sum((a-tm)*(b-Lm) for a,b in zip(t,L))/sum((a-tm)**2 for a in t)   # hbar/s
    r['slope']=slope*1e-9; r['dLmax']=max(abs(x-L[0]) for x in L)
    r['tcross']=0.5/abs(slope)*1e9 if slope!=0 else float('inf')
    r['wmin']=min(r['w']); r['wmean']=st.mean(r['w']); r['wmax']=max(r['w'])
    P(f"{r['i']:3d} {r['ch']} {len(t):4d} {t[-1]*1e12:9.1f} {L[0]:.5f} {r['dLmax']:.2e}  {r['slope']:+.3e}      {r['tcross']:.3e}     [{min(r['n']):.6f},{max(r['n']):.6f}] {r['wmin']:.4f}  {r['wmean']:.4f}  {r['wmax']:.4f}")
P("\nSUMMARY per channel")
for ch in ['p-Ps','o-Ps']:
    R=[r for r in res if r['ch']==ch and 'slope' in r]
    if not R: continue
    P(f"{ch}: N={len(R)}  max|dL| over all {max(r['dLmax'] for r in R):.2e} hbar;  t(0.5 hbar) median {st.median(r['tcross'] for r in R):.3e} ns, min {min(r['tcross'] for r in R):.3e} ns;"
      f" > 142 ns: {sum(r['tcross']>142 for r in R)}/{len(R)};  w min<0.8: {sum(r['wmin']<0.8 for r in R)}/{len(R)};  max w {max(r['wmax'] for r in R):.3e};  median time-mean w {st.median(r['wmean'] for r in R):.4f}")
