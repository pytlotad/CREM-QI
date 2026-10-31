# Audit 665: the engine's l >= 1 sum rule (e_all) with other orbital angular momentum lengths, offline, hydrogen (infinite nuclear mass; every variant
# scales as mu, so the factor is the same for p-Ps and o-Ps).  usage: python3 rate665.py   (run from this directory)
import numpy as np
src=open('../section-560/tempo_e1_punkt_srodkowy.py').read().split('if __name__ == "__main__":')[0]
g={}; exec(src,g)
rate_qm,rate_rule,channels,C,AU=g['rate_qm'],g['rate_rule'],g['channels'],g['C'],g['AU_RATE']
def e_all(n,L,N=1<<20):          # = 585 / 664 e_all with the orbit's L given directly (units hbar)
    a=n*n; e=np.sqrt(max(0.0,1-(L/n)**2)); b=a*np.sqrt(1-e*e)
    M=(np.arange(N))*2*np.pi/N; u=M.copy() if e<1e-12 else np.full_like(M,np.pi)
    for _ in range(200): u-= (u-e*np.sin(u)-M)/(1-e*np.cos(u))
    z=a*(np.cos(u)-e)+1j*b*np.sin(u); c=np.fft.fft(z)/N
    k=np.arange(1,N//2); Om=1.0/n**3
    return np.sum(2/3*np.abs(c[k])**2*(k*Om)**3)/C**3
xg,wg=np.polynomial.legendre.leggauss(8)
ns=lambda A: 1e9/(A*AU)
names="spdfghiklm"; R={c:[] for c in ('S0','S1','S2')}
print('state  tau_QM[ns]   tau/QM:  S0      S1      S2')
for n in range(2,11):
    for l in range(1,n):
        q=sum(rate_qm(n,l,n2,l2) for n2,l2 in channels(n,l))
        ru=sum(rate_rule(n,l,n2,l2,'RBw') for n2,l2 in channels(n,l) if l2==l+1)
        A={'S0':e_all(n,l+0.5)+ru,'S1':e_all(n,np.sqrt(l*(l+1)))+ru,
           'S2':sum(0.5*w*e_all(n,l+0.5+0.5*x) for x,w in zip(xg,wg))+ru}
        row=f'{n}{names[l]:<4} {ns(q):10.3f}        '
        for cn in R: t=q/A[cn]; R[cn].append(t); row+=f'{t:7.4f} '
        print(row,flush=True)
print()
s0=np.array(R['S0'])
for cn,v in R.items():
    v=np.array(v); mx=np.max(np.abs(v-1))
    nw=np.sum(np.abs(np.log(v))<=np.abs(np.log(s0))+1e-12)
    print(f'{cn}: max |tau/QM-1| {mx:.4f}, 2p {v[0]:.4f}; H665a {mx<=0.01 and abs(v[0]-1)<=0.01}; not worse than S0 in {nw}/{len(v)}, '
          f'H665b {cn!="S0" and nw>=0.9*len(v) and mx<np.max(np.abs(s0-1))}')
