# Audit 506: emission time of the 2p -> 1s photon vs spin orientation, per variant.  python3 ana506.py <dir>
import math, sys, glob, os
def V(a): return [float(x) for x in a[:3]]
def dot(a,b): return sum(x*y for x,y in zip(a,b))
def nrm(a): return math.sqrt(dot(a,a))
def corr(x,y):
    mx=sum(x)/len(x); my=sum(y)/len(y); sx=math.sqrt(sum((a-mx)**2 for a in x)); sy=math.sqrt(sum((b-my)**2 for b in y))
    return sum((a-mx)*(b-my) for a,b in zip(x,y))/(sx*sy) if sx>0 and sy>0 else float('nan')
d=sys.argv[1]
for v in ('default','NO_DIPOLE_FORCE','NO_THOMAS_BACKREACTION','NO_HIDDEN_MOMENTUM_FORCE','NO_SPIN_ENERGY_EXCHANGE'):
    rows=[]; idx=None; ch=None
    for f in sorted(glob.glob(os.path.join(d,v+'.*.err'))):
        for l in open(f):
            if l.startswith('START'): idx=int(l.split()[2])
            elif l.startswith('BORN'): ch=l.split()[-1]
            elif l.startswith('CREM_LEDGER '):
                t=l.split(); b=[k for k,x in enumerate(t) if x=='|']; s=t[b[0]+2:]
                m1=V(s[s.index('m1')+1:]); m2=V(s[s.index('m2')+1:]); ax=V(s[s.index('ax')+1:])
                rows.append(dict(i=idx,t=float(t[2]),ceL=-dot(m2,ax)/(nrm(m2)*nrm(ax)),cpL=dot(m1,ax)/(nrm(m1)*nrm(ax)),c12=-dot(m1,m2)/(nrm(m1)*nrm(m2)),ch=ch))
    if not rows: print(f'{v}: no data'); continue
    ts=[r['t'] for r in rows]; m=sum(ts)/len(ts)
    print(f"{v:26s} n {len(rows):2d}: t mean {m*1e9:.6f} ns, spread (max-min)/mean {(max(ts)-min(ts))/m:.2e}; corr(t, cos(s_e,L)) {corr(ts,[r['ceL'] for r in rows]):+.3f}, "
          f"corr(t, cos(s_p,L)) {corr(ts,[r['cpL'] for r in rows]):+.3f}, corr(t, cos(s1,s2)) {corr(ts,[r['c12'] for r in rows]):+.3f}")
    if v=='default':
        rs=sorted(rows,key=lambda r:r['ceL']); print('   extremes (idx, cos(s_e,L), t ns, channel):',[(r['i'],round(r['ceL'],3),round(r['t']*1e9,6),r['ch']) for r in (rs[0],rs[1],rs[-2],rs[-1])])
        # slope dt/dcos
        x=[r['ceL'] for r in rows]; mx=sum(x)/len(x); slope=sum((a-mx)*(b-m) for a,b in zip(x,ts))/sum((a-mx)**2 for a in x)
        print('   slope dt/dcos(s_e,L) = %.3e s (relative %.2e)'%(slope,slope/m))
