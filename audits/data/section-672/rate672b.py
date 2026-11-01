# Audit 672, variant V2: the orbit on which both balances hold, P^+/count = E_bar (the engine's mean photon energy), offline, hydrogen.
# usage: python3 rate672b.py   (run from this directory; uses rate672.py's functions)
import numpy as np
src=open('rate672.py').read().split('ns=lambda')[0]; g={}; exec(src,g)
harmonics,ebar,rate_qm,rate_rule,channels,AU=g['harmonics'],g['ebar'],g['rate_qm'],g['rate_rule'],g['channels'],g['AU']
ns=lambda A: 1e9/(A*AU)
def ratio(n,L,N):
    cnt,pw,_,_=harmonics(n,L,N); return pw/cnt, cnt
names="spdfghiklm"; T0=[]; T2=[]
print('state  tau_QM[ns]  E_bar/Omega  <kOmega>/E_bar(Langer)  L*/(l+1/2)  tau/QM: C0      V2')
for n in range(2,11):
    for l in range(1,n):
        q=sum(rate_qm(n,l,n2,l2) for n2,l2 in channels(n,l))
        ru=sum(rate_rule(n,l,n2,l2,'RBw') for n2,l2 in channels(n,l) if l2==l+1)
        Eb,_=ebar(n,l); Om=1.0/n**3
        r0,c0=ratio(n,l+0.5,1<<20)
        f=lambda L: ratio(n,L,1<<16)[0]-Eb
        lo,hi=0.02,float(n); flo,fhi=f(lo),f(hi)
        if flo*fhi>0: Ls=None
        else:
            for _ in range(40):
                mid=0.5*(lo+hi); fm=f(mid)
                if fm*flo>0: lo,flo=mid,fm
                else: hi=mid
            Ls=0.5*(lo+hi)
        cs=ratio(n,Ls,1<<20)[1] if Ls else float('nan')
        t0=q/(c0+ru); t2=q/(cs+ru); T0.append(t0); T2.append(t2)
        print(f'{n}{names[l]:<4} {ns(q):10.3f}   {Eb/Om:9.4f}   {r0/Eb:12.4f}            '+(f'{Ls/(l+0.5):8.4f}' if Ls else '    none')+f'      {t0:7.4f} {t2:7.4f}',flush=True)
T0=np.array(T0); T2=np.array(T2); ok=np.isfinite(T2)
mx=np.max(np.abs(T2[ok]-1)) if ok.any() else float('nan'); nw=np.sum(np.abs(np.log(T2[ok]))<=np.abs(np.log(T0[ok]))+1e-12)
print(f'\nV2: solutions {ok.sum()}/{len(T2)}; max |tau/QM-1| {mx:.4f}, 2p {T2[0]:.4f}; H672c {bool(ok.all() and mx<=0.01 and abs(T2[0]-1)<=0.01)};'
      f' not worse than C0 in {nw}/{len(T2)}, H672d {bool(ok.all() and nw>=0.9*len(T2) and mx<np.max(np.abs(T0-1)))}')
