# Audit 598: per-mode vs net radiation torque; reversal <=> parity; E1/E2 channel criteria; pair and quadrupole torques.
import numpy as np
def kepler(e,N):
    M=np.arange(N)*2*np.pi/N; u=np.full_like(M,np.pi)
    for _ in range(200): u-=(u-e*np.sin(u)-M)/(1-e*np.cos(u))
    return M,u
print("T1  per-mode torque (Kepler, a = 1, mean motion 1): sum_s P_{s,m}/s for m = +1 (co), -1 (counter)")
for e in (0.9,0.99,0.999):
    N=1<<20; M,u=kepler(e,N); z=(np.cos(u)-e)+1j*np.sqrt(1-e*e)*np.sin(u)
    c=np.fft.fft(z)/N; k=np.arange(1,N//2)
    co=np.sum(np.abs(c[k])**2*k**4/k); cn=np.sum(np.abs(c[N-k])**2*k**4/k)   # P_s ~ s^4 |c|^2, torque ~ P/s
    print(f"    e={e}: counter/co = {cn/co:.4f};  net/co = {(co-cn)/co:.4f}")
bad=0
for l in range(0,41):
    for m,par in ((1,1),(-1,1),(2,0),(-2,0)):      # par: required parity of Delta l (1 odd for E1, 0 even for E2)
        L=l+0.5; Lp=abs(L-m); lp=round(Lp-0.5); dl=lp-l
        rev=m>L; wrong=(abs(dl)%2)!=par
        bad+=rev!=wrong
print(f"T2  reversal <=> wrong parity: {bad} mismatches (l = 0..40, m = +-1, +-2)")
def qmE2(l,lp): return lp>=0 and lp in (l-2,l,l+2) and not (l==0 and lp==0)
def qmE1(l,lp): return lp>=0 and lp in (l-1,l+1)
def tri(l,lp,j): L,Lp,Lph=l+0.5,lp+0.5,j+0.5; return abs(L-Lph)<=Lp<=L+Lph
def engE2(l,m):
    t=(2.0 if l>=2 else -2.0) if m==2 else ((-2.0 if l==0 else 0.0) if m==0 else -2.0)
    return round(abs(l+0.5-t)-0.5)
res={'a':0,'b':0,'c':0,'d':0,'cE1':0}; extra={'a':[],'b':[]}
for n in range(1,21):
  for l in range(n):
    for n2 in range(1,n+1):
      for lp in range(0,n2):
        q2=qmE2(l,lp); q1=qmE1(l,lp)
        va=set(); 
        for m in (2,0,-2):
            L=l+0.5
            if m>L: continue
            va.add(round(abs(L-m)-0.5))
        a=lp in va; b=a and (abs(lp-l)%2==0); c=tri(l,lp,2) and (abs(lp-l)%2==0); d=lp in {engE2(l,m) for m in (2,0,-2)}
        for key,v in (('a',a),('b',b),('c',c),('d',d)):
            if v!=q2:
                res[key]+=1
                if key in extra and len(extra[key])<4: extra[key].append((l,lp,v,q2))
        c1=tri(l,lp,1) and (abs(lp-l)%2==1)
        res['cE1']+=c1!=q1
print("T3  E2 channel sets vs QM (mismatches): (a) vector + no reversal", res['a'], extra['a'], "; (b) + parity", res['b'], extra['b'],
      "; (c) triangle with rule-Q lengths + parity", res['c'], "; (d) engine rule", res['d'], ";  E1 (c):", res['cE1'])
# T4 (i) pair relative dipole d = q r: r x d''' = q r x a' = -q L/r^3 (same identity as 597), (ii) quadrupole torque sign
for e in (0.1,0.5,0.9,0.99,0.999):
    N=1<<16; M,u=kepler(e,N); x=np.cos(u)-e; y=np.sqrt(1-e*e)*np.sin(u)
    def deriv(f,n):   # spectral derivative in time (mean motion 1)
        F=np.fft.fft(f); w=1j*np.fft.fftfreq(N,1.0/N); return np.real(np.fft.ifft(F*w**n))
    Qxx=x*x-(x*x+y*y)/3; Qyy=y*y-(x*x+y*y)/3; Qxy=x*y
    t2=lambda q: deriv(q,2); t3=lambda q: deriv(q,3)
    # torque_z ~ -eps_zjk Q''_jl Q'''_kl = -(Q''_xl Q'''_yl - Q''_yl Q'''_xl), l in x, y (z-row of Q is diagonal only)
    tz=-np.mean(t2(Qxx)*t3(Qxy)+t2(Qxy)*t3(Qyy)-t2(Qxy)*t3(Qxx)-t2(Qyy)*t3(Qxy))
    print(f"T4  e={e}: quadrupole torque along L^ (L_z > 0): {tz:+.3e}  (negative = L decreases)")
