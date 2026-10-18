# Audit 504: Zeeman photon (import 20) -- 2p -> 1s in 10 T vs B = 0 (same seeds, section-503/B0.err).  python3 ana504.py Z.err B0.err
import math, sys
hbar=1.054571817e-34; eV=1.602176634e-19; muB=9.2740100783e-24; mp=1.67262192369e-27; me=9.1093837015e-31
def V(a): return [float(x) for x in a[:3]]
def dot(a,b): return sum(x*y for x,y in zip(a,b))
def nrm(a): return math.sqrt(dot(a,a))
def unit(a): n=nrm(a); return [x/n for x in a]
def load(f):
    out={}; idx=None; B=None
    for l in open(f):
        if l.startswith('BFIELD'): B=[float(x) for x in l.split()[3:6]]
        elif l.startswith('START'): idx=int(l.split()[2])
        elif l.startswith('CREM_LEDGER ') and idx is not None and idx not in out:
            t=l.split(); b=[k for k,x in enumerate(t) if x=='|']; sb=t[b[0]+2:]; sa=t[b[1]+2:]; g=lambda s,k: s[s.index(k)+1:]
            out[idx]=dict(B=B,t=float(t[2]),Eg=float(t[t.index('Eg')+1]),dir=V(t[t.index('dir')+1:]),hel=float(t[t.index('hel')+1]),
                sph=V(t[t.index('sphoton')+1:]),LB=float(g(sb,'L')[0]),axB=V(g(sb,'ax')),LA=float(g(sa,'L')[0]),axA=V(g(sa,'ax')),
                UF=float(g(sb,'UF')[0])-float(g(sa,'UF')[0]))
    return out
Z=load(sys.argv[1]); O=load(sys.argv[2]); keys=sorted(set(Z)&set(O))
wL=muB*mp/(mp+me)   # hbar omega_L per tesla (reduced-mass Larmor), J/T
rows=[]
for k in keys:
    z,o=Z[k],O[k]; Bm=nrm(z['B']); bh=unit(z['B'])
    m=dot(z['sph'],bh)/hbar; c=dot(unit(z['axB']),bh)
    resid=(z['Eg']-o['Eg']-round(m)*wL*Bm)/(muB*Bm)
    rows.append(dict(k=k,m=m,c=c,dE=(z['Eg']-o['Eg'])/(muB*Bm),resid=resid,ck=dot(unit(z['dir']),bh),hel=z['hel'],t=z['t'],
                     LAz=dot(unit(z['axA']),bh)*z['LA']/hbar,LBz=c*z['LB']/hbar))
n=len(rows); print(f'paired photons: {n}')
print('P1 max |m_gamma - round(m_gamma)| = %.2e; populations -1/0/+1: %s'%(max(abs(r['m']-round(r['m'])) for r in rows),
      [sum(1 for r in rows if round(r['m'])==v) for v in (-1,0,1)]))
cnt=[sum(1 for r in rows if round(r['m'])==v) for v in (-1,0,1)]; chi=sum((x-n/3)**2/(n/3) for x in cnt)
print('   chi^2 vs 1/3 each = %.2f (2 dof; p = %.3f); basin check: Delta m = sign rule of L.B^/hbar (cos > 1/3 -> +1 ...): %d/%d'%(chi,math.exp(-chi/2),
      sum(1 for r in rows if round(r['m'])==(1 if r['c']>1/3 else (-1 if r['c']<-1/3 else 0))),n))
print('   L.B^/hbar before -> after: '+' '.join(f"{r['LBz']:+.2f}->{r['LAz']:+.2f}" for r in rows[:8]))
print('P2 energy shift / (mu_B B) by Delta m:')
for v in (-1,0,1):
    x=[r['dE'] for r in rows if round(r['m'])==v]
    if x: print(f'   Delta m {v:+d}: n {len(x):2d}, shift {min(x):+.4f} .. {max(x):+.4f} (expected {v*mp/(mp+me):+.5f} + spin-level part)')
print('   residual after removing Delta m hbar omega_L: max |.| = %.4f mu_B B'%max(abs(r['resid']) for r in rows))
for v,lab,exp in ((0,'pi',0.2),(1,'sigma',0.4)):
    x=[r['ck']**2 for r in rows if (round(r['m'])==0)==(v==0)]
    if x: print(f'P3 {lab}: n {len(x)}, <cos^2(k,B)> = {sum(x)/len(x):.3f} (expected {exp}; sd of mean ~ {math.sqrt((0.0571 if v==0 else 0.1086)/len(x)):.3f})')
s=[r for r in rows if round(r['m'])!=0]
print('P4 sigma photons: helicity sign = sign(Delta m cos(k,B)) in %d/%d (expected ~ (1+|c|)^2/(2(1+c^2)) on average > 1/2)'%(
      sum(1 for r in s if r['hel']*round(r['m'])*r['ck']>0),len(s)))
print('P5 emission time %.4f - %.4f ns'%(min(r['t'] for r in rows)*1e9,max(r['t'] for r in rows)*1e9))
