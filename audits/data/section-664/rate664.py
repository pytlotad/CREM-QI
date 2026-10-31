# Audit 664: E1 rate rules for l >= 1 against exact QM, offline, hydrogen (infinite nuclear mass; the rate of every rule scales as mu, so the
# factor is the same for p-Ps and o-Ps).  usage: python3 rate664.py   (run from this directory)
import numpy as np
src=open('../section-560/tempo_e1_punkt_srodkowy.py').read().split('if __name__ == "__main__":')[0]
g={}; exec(src,g)
rate_qm,rate_rule,channels,fc,omega,C,AU=g['rate_qm'],g['rate_rule'],g['channels'],g['fourier_components'],g['omega'],g['C'],g['AU_RATE']
def e_all(n,l,N=1<<20,parts=False):          # = 585 e_all (the engine's l >= 1 count)
    a=n*n; L=l+0.5; e=np.sqrt(max(0.0,1-(L/n)**2)); b=a*np.sqrt(1-e*e)
    M=(np.arange(N))*2*np.pi/N; u=np.full_like(M,np.pi)
    for _ in range(200): u-= (u-e*np.sin(u)-M)/(1-e*np.cos(u))
    z=a*(np.cos(u)-e)+1j*b*np.sin(u); c=np.fft.fft(z)/N
    k=np.arange(1,N//2); Om=1.0/n**3; t=2/3*np.abs(c[k])**2*(k*Om)**3/C**3
    return (t[0],np.sum(t[1:])) if parts else np.sum(t)
def mode(n_orb,L,s,down): cp,cm=fc(n_orb,L,s,N=16384); return cp if down else cm
def cand(name,n,l,n2,l2):
    s=n-n2; down=l2==l-1
    c=mode(n,l+0.5,s,down)
    if name in ('C3','C4'): c=np.sqrt(c*mode(n2,l2+0.5,s,down))
    w=omega(n,n2) if name in ('C2','C4') else s/n**3
    return 2/3*c**2*w**3/C**3
ns=lambda A: 1e9/(A*AU)
# step A
k1,k2=e_all(2,1,parts=True); qm2=rate_qm(2,1,1,0)
print(f'STEP A, 2p (e^2 = 0.4375): E_all = {k1+k2:.6e} au, harmonic k=1 {k1:.6e} ({100*k1/(k1+k2):.2f} %), k>=2 {k2:.6e} ({100*k2/(k1+k2):.2f} %); '
      f'QM {qm2:.6e}; E_all/QM {(k1+k2)/qm2:.4f}, k=1 only / QM {k1/qm2:.4f}')
print(f'  tau 2p: QM {ns(qm2):.4f} ns, E_all {ns(k1+k2):.4f} ns (engine 663: 1.6748 ns)')
names="spdfghiklm"; R={c:[] for c in ('C0','C1','C2','C3','C4')}
print('\nstate  tau_QM[ns]   tau/QM:  C0      C1      C2      C3      C4')
for n in range(2,11):
    for l in range(1,n):
        q=sum(rate_qm(n,l,n2,l2) for n2,l2 in channels(n,l))
        ru=sum(rate_rule(n,l,n2,l2,'RBw') for n2,l2 in channels(n,l) if l2==l+1)
        A={'C0':e_all(n,l)+ru}
        for cn in ('C1','C2','C3','C4'): A[cn]=sum(cand(cn,n,l,n2,l2) for n2,l2 in channels(n,l))
        row=f'{n}{names[l]:<4} {ns(q):10.3f}        '
        for cn in R: t=q/A[cn]; R[cn].append(t); row+=f'{t:7.3f} '
        print(row)
print()
c0=np.array(R['C0'])
for cn,v in R.items():
    v=np.array(v); mx=np.max(np.abs(v-1)); p2=abs(v[0]-1)
    nw=np.sum(np.abs(np.log(v))<=np.abs(np.log(c0))+1e-12)
    print(f'{cn}: max |tau/QM-1| {mx:.4f}, 2p {v[0]:.4f}; H664a {mx<=0.01 and p2<=0.01}; not worse than C0 in {nw}/{len(v)}, '
          f'H664b {cn!="C0" and nw>=0.9*len(v) and mx<np.max(np.abs(c0-1))}')
