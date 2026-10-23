# Audit 600: Kramers path average of the classical Fourier amplitude along the transition (hydrogen, atomic units).
import math, numpy as np
src=open('../section-560/tempo_e1_punkt_srodkowy.py').read().split('if __name__ == "__main__":')[0]
g={}; exec(src,g)
fc,rate_qm,channels,C=g['fourier_components'],g['rate_qm'],g['channels'],g['C']
xg,wg=np.polynomial.legendre.leggauss(16); tg=0.5*(xg+1); wg=0.5*wg
cache={}
def amps(n,l,n2,l2):
    key=(n,l,n2,l2)
    if key in cache: return cache[key]
    s=n-n2; out=[]
    for t in tg:
        nt=n-t*s; Lt=(l+0.5)+t*(l2-l)
        cp,cm=fc(nt,Lt,s,N=8192); out.append((cp if l2==l-1 else cm, s/nt**3))
    cache[key]=out; return out
def rate(v,n,l,n2,l2):
    A=amps(n,l,n2,l2); w=2*max(l,l2)/(2*l+1); wx=0.5*(1/n2**2-1/n**2)
    c=np.array([a[0] for a in A]); om=np.array([a[1] for a in A])
    if v=='K1': return w*2/3*(wg@c)**2*wx**3/C**3
    if v=='K2': return w*2/3*(wg@(c*c))*wx**3/C**3
    if v=='K3': return w*2/3*(wg@(c*om**1.5))**2/C**3
for v in ('K1','K2','K3'):
    lr=[math.log(rate(v,n,l,n2,l2)/rate_qm(n,l,n2,l2)) for n in range(2,11) for l in range(n) for n2,l2 in channels(n,l)]
    big=[rate(v,n,1,1,0)/rate_qm(n,1,1,0) for n in range(3,9)]+[rate(v,n,0,2,1)/rate_qm(n,0,2,1) for n in range(3,11)]
    lt=[]
    for n in range(2,9):
        for l in range(min(n,4)):
            if (n,l)==(2,0): continue
            q=sum(rate_qm(n,l,n2,l2) for n2,l2 in channels(n,l)); r=sum(rate(v,n,l,n2,l2) for n2,l2 in channels(n,l))
            lt.append((q/r,f"{n}{'spdf'[l]}"))
    worst=max(lt,key=lambda t:abs(t[0]-1))
    print(f"{v}: channels rms|ln| {math.sqrt(sum(x*x for x in lr)/len(lr)):.3f} max|ln| {max(abs(x) for x in lr):.3f};  largest-Delta-n {min(big):.3f}-{max(big):.3f};  lifetimes worst {worst[1]} {worst[0]:.3f}",flush=True)
    print("    "+"  ".join(f"{name} {x:.3f}" for x,name in lt))
