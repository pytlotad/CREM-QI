# Audit 612: the engine's np Delta l = -1 rate decomposed by harmonic/channel vs exact QM (hydrogen).
import numpy as np
src=open('../section-560/tempo_e1_punkt_srodkowy.py').read().split('if __name__ == "__main__":')[0]
g={}; exec(src,g)
rate_qm,rate_rule,C=g['rate_qm'],g['rate_rule'],g['C']
def harmonics(n,l,N=1<<20):
    a=n*n; L=l+0.5; e=np.sqrt(max(0.0,1-(L/n)**2)); b=a*np.sqrt(1-e*e)
    M=np.arange(N)*2*np.pi/N; u=np.full_like(M,np.pi)
    for _ in range(200): u-=(u-e*np.sin(u)-M)/(1-e*np.cos(u))
    z=a*(np.cos(u)-e)+1j*b*np.sin(u); c=np.fft.fft(z)/N; k=np.arange(1,N//2); Om=1.0/n**3
    return k, 2/3*np.abs(c[k])**2*(k*Om)**3/C**3            # photon-count rate of harmonic k (co-rotating)
print("n   tau_eng/QM   Gamma deficit by channel (A_eng - A_QM)/Gamma_QM: n' = n-1 ... 1 (1 = all s >= n-1), up (n'd)   | beyond-grid share (s >= n)")
for n in range(2,11):
    l=1; k,r=harmonics(n,l)
    down_qm={n2:rate_qm(n,l,n2,0) for n2 in range(1,n)}
    up_qm={n2:rate_qm(n,l,n2,2) for n2 in range(3,n)}; up_eng={n2:rate_rule(n,l,n2,2,'RBw') for n2 in up_qm}
    eng={}
    for n2 in range(2,n): eng[n2]=r[k==n-n2].sum()
    eng[1]=r[k>=n-1].sum()
    G=sum(down_qm.values())+sum(up_qm.values()); Ge=sum(eng.values())+sum(up_eng.values())
    parts=[(n2,(eng[n2]-down_qm[n2])/G) for n2 in range(n-1,0,-1)]
    upd=(sum(up_eng.values())-sum(up_qm.values()))/G
    beyond=r[k>=n].sum()/Ge
    print(f"{n:2d}   {G/Ge:.4f}     "+" ".join(f"{n2}s:{d:+.4f}" for n2,d in parts)+f"  up:{upd:+.4f}   | {beyond:.3f}")
print("\nper-channel ratio A_eng/A_QM for np -> n's (n' = n-1 ... 1):")
for n in range(2,11):
    k,r=harmonics(n,1)
    row=[]
    for n2 in range(n-1,0,-1):
        ae=r[k==n-n2].sum() if n2>1 else r[k>=n-1].sum(); row.append(f"{n2}s {ae/rate_qm(n,1,n2,0):.3f}")
    print(f"  {n}p: "+"  ".join(row)+f"   (1s from s = n-1 only: {r[k==n-1].sum()/rate_qm(n,1,1,0):.3f})")
