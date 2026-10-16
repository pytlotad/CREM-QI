# Audit 492 B/C: the engine's spin vectors in the author's family sigma_1,2 = a S^ +- u (sigma = 1 for the engine's hbar/2).
# B: 492 (default: free draw + Born channel), C: 490 (CREM_NO_BORN_SPIN_CHANNEL, free spins).  At every photon (after side):
# a = |s1 + s2|/2 = |S|/hbar, |u| = |s1 - s2|/2, cos12, and the engine vectors' own orbit-averaged tensor on the planar orbit
# (isotropic r^-3 weight, 491a): <S12> = 1.5 (s1_perp . s2_perp) - s1 . s2, perp = in the orbit plane.
# Channel = sign of s1.s2 (the engine's spinsAntiparallel); 490 has no Born channel, so there it is only the sign.
# python3 ana492.py <492 Ps.err files...>
import math, sys, glob
hbar=1.054571817e-34
def V(a): return [float(x) for x in a[:3]]
def dot(a,b): return sum(x*y for x,y in zip(a,b))
def norm(a): return math.sqrt(dot(a,a))
def unit(a): n=norm(a); return [x/n for x in a]
def add(a,b): return [x+y for x,y in zip(a,b)]
def mul(a,k): return [x*k for x in a]
def perp(v,n): return add(v,mul(n,-dot(v,n)))
def rows(files):
    out=[]
    for f in files:
        for line in open(f):
            if not line.startswith('CREM_LEDGER '): continue
            t=line.split(); seg=t[[k for k,x in enumerate(t) if x=='|'][1]+2:]
            g=lambda key: seg.index(key)
            m1=V(seg[g('m1')+1:]); m2=V(seg[g('m2')+1:])
            # channel as the engine reads it (spinsAntiparallel): the ledger's 'singlet' field is gSingletChannelEligible, not the channel
            out.append(dict(singlet=int(dot(mul(unit(m1),-1),unit(m2))<0.0),L=float(seg[g('L')+1])/hbar,ax=unit(V(seg[g('ax')+1:])),m1=V(seg[g('m1')+1:]),m2=V(seg[g('m2')+1:])))
    return out
def report(name,files):
    R=rows(files); sg=(-1,1)   # e- moment antiparallel to its spin, e+ parallel
    print(f'== {name}: photons {len(R)} (triplet {sum(1 for r in R if not r["singlet"])}, singlet {sum(r["singlet"] for r in R)})')
    for ch,lab in ((0,'triplet'),(1,'singlet')):
        X=[r for r in R if r['singlet']==ch]
        if not X: continue
        dev=0; A=[]; U=[]; S12=[]; S12P=[]
        for r in X:
            s1=mul(unit(r['m1']),sg[0]); s2=mul(unit(r['m2']),sg[1]); c=dot(s1,s2)
            dev=max(dev,abs(1-c) if ch==0 else abs(1+c))
            A.append(norm(add(s1,s2))/2); U.append(norm(add(s1,mul(s2,-1)))/2)
            v=1.5*dot(perp(s1,r['ax']),perp(s2,r['ax']))-c; S12.append(v)
            if r['L']>0.99: S12P.append(v)
        print(f'   {lab}: max |1 -+ cos12| {dev:.2e};  a = |S|/hbar in [{min(A):.3f}, {max(A):.3f}];  |u| in [{min(U):.3f}, {max(U):.3f}]')
        print(f'      engine-vector <S12> all photons [{min(S12):+.3f}, {max(S12):+.3f}];  after-states L >= hbar: {len(S12P)}'
              +(f', [{min(S12P):+.3f}, {max(S12P):+.3f}], mean {sum(S12P)/len(S12P):+.3f}' if S12P else ''))
report('492 Ps default (Born channel)',sys.argv[1:])
report('490 Ps free spins (no Born channel)',['../section-490/Ps.err'])
