# Audit 672: E1 count from the energy balance (classical co-rotating power / the engine's mean photon energy) against exact QM, offline, hydrogen.
# usage: python3 rate672.py   (run from this directory)
import numpy as np
src=open('../section-560/tempo_e1_punkt_srodkowy.py').read().split('if __name__ == "__main__":')[0]
g={}; exec(src,g)
rate_qm,rate_rule,channels,fc,omega,C,AU=g['rate_qm'],g['rate_rule'],g['channels'],g['fourier_components'],g['omega'],g['C'],g['AU_RATE']
def harmonics(n,L,N=1<<20):
    a=n*n; e=np.sqrt(max(0.0,1-(L/n)**2)); b=a*np.sqrt(1-e*e)
    M=(np.arange(N))*2*np.pi/N; u=M.copy() if e<1e-12 else np.full_like(M,np.pi)
    for _ in range(200): u-= (u-e*np.sin(u)-M)/(1-e*np.cos(u))
    z=a*(np.cos(u)-e)+1j*b*np.sin(u); c=np.fft.fft(z)/N
    k=np.arange(1,N//2); Om=1.0/n**3; base=2/3*np.abs(c[k])**2/C**3
    return np.sum(base*(k*Om)**3), np.sum(base*(k*Om)**4), k, base*(k*Om)**3   # count (L balance), power, harmonic count terms
def ebar(n,l):
    w=[]; gaps=[]
    for n2 in range(l,n):
        s=n-n2; nm=np.sqrt(n*n2); cp,_=fc(nm,l,s,N=16384); w.append(cp**2*(s/nm**3)**3); gaps.append(omega(n,n2))
    w=np.array(w)/sum(w); return float(w@np.array(gaps)), w
ns=lambda A: 1e9/(A*AU)
cnt,pw,k,terms=harmonics(2,1.5); Eb,_=ebar(2,1)
keff3=np.sum(terms*k)/np.sum(terms)
print(f'STEP A 2p: L-balance count {cnt:.6e}, power/E_bar {pw/Eb:.6e} au; E_bar = {Eb:.6f} au (= gap 3/8), Omega = {1/8:.6f}; '
      f'V1/C0 = Omega <k> / E_bar = {pw/Eb/cnt:.4f} (count-weighted <k> = {keff3:.4f}); QM {rate_qm(2,1,1,0):.6e}')
names="spdfghiklm"; R={'C0':[],'V1':[]}
print('\nstate  tau_QM[ns]   tau/QM:  C0      V1')
for n in range(2,11):
    for l in range(1,n):
        q=sum(rate_qm(n,l,n2,l2) for n2,l2 in channels(n,l))
        ru=sum(rate_rule(n,l,n2,l2,'RBw') for n2,l2 in channels(n,l) if l2==l+1)
        cnt,pw,_,_=harmonics(n,l+0.5); Eb,_=ebar(n,l)
        A={'C0':cnt+ru,'V1':pw/Eb+ru}
        row=f'{n}{names[l]:<4} {ns(q):10.3f}        '
        for cn in R: t=q/A[cn]; R[cn].append(t); row+=f'{t:7.4f} '
        print(row,flush=True)
print()
c0=np.array(R['C0'])
for cn,v in R.items():
    v=np.array(v); mx=np.max(np.abs(v-1)); nw=np.sum(np.abs(np.log(v))<=np.abs(np.log(c0))+1e-12)
    print(f'{cn}: max |tau/QM-1| {mx:.4f}, 2p {v[0]:.4f}; H672a {mx<=0.01 and abs(v[0]-1)<=0.01}; not worse than C0 in {nw}/{len(v)}, '
          f'H672b {cn!="C0" and nw>=0.9*len(v) and mx<np.max(np.abs(c0-1))}')
