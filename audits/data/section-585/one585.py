# Audit 585, addendum 2: can one place serve both directions of Delta l?
import numpy as np
src=open('../section-560/tempo_e1_punkt_srodkowy.py').read().split('if __name__ == "__main__":')[0]
g={}; exec(src,g)
rate_qm, rate_rule, channels, C = g['rate_qm'], g['rate_rule'], g['channels'], g['C']
fc=g['fourier_components']
def e_all(n,l,m,N=1<<20):
    a=n*n; L=l+0.5; e=np.sqrt(max(0.0,1-(L/n)**2)); b=a*np.sqrt(1-e*e)
    M=np.arange(N)*2*np.pi/N; u=np.full_like(M,np.pi)
    for _ in range(200): u-=(u-e*np.sin(u)-M)/(1-e*np.cos(u))
    c=np.fft.fft(z:=a*(np.cos(u)-e)+1j*b*np.sin(u))/N
    k=np.arange(1,N//2); Om=1.0/n**3
    ck=c[k] if m>0 else c[N-k]
    return np.sum(2/3*np.abs(ck)**2*(k*Om)**3)/C**3
names="spdfghijk"
print("(i) Delta l = -1 at the midpoint with weight 2l/(2l+1), per channel / QM, then lifetime of the down channels alone / QM")
for n,l in [(3,1),(5,1),(8,1),(5,2),(8,2),(8,4),(10,6)]:
    cells=[f"{n2}{names[l-1]} {rate_rule(n,l,n2,l-1,'RBw')/rate_qm(n,l,n2,l-1):.2f}" for n2 in range(n-1,l-1,-1)]
    qd=sum(rate_qm(n,l,n2,l-1) for n2 in range(l,n)); rd=sum(rate_rule(n,l,n2,l-1,'RBw') for n2 in range(l,n))
    print(f"  {n}{names[l]}: tau_down {qd/rd:.3f} | "+' '.join(cells))
print("\n(ii) Delta l = +1 at the initial Langer orbit (counter-rotating), up-share / QM:  per channel w=1 | per channel w=2l_>/(2l+1) | all harmonics (E_all m=-1)  [midpoint RBw]")
for n,l in [(5,1),(7,1),(10,1),(7,2),(10,2),(10,4),(10,6)]:
    qd=sum(rate_qm(n,l,n2,l-1) for n2 in range(l,n)); qu=sum(rate_qm(n,l,n2,l+1) for n2 in range(l+2,n))
    AE=e_all(n,l,+1)
    def up(f): return f/(AE+f)/(qu/(qd+qu))
    r0=0.0
    for n2 in range(l+2,n):
        s=n-n2; cm=fc(n,l+0.5,s)[1]; r0+=2/3*cm**2*(s/n**3)**3/C**3
    w=2*(l+1)/(2*l+1); allm=e_all(n,l,-1)
    rb=sum(rate_rule(n,l,n2,l+1,'RBw') for n2 in range(l+2,n))
    print(f"  {n}{names[l]}: {up(r0):.3f} | {up(w*r0):.3f} | {up(allm):.3f}  [{up(rb):.3f}]")
