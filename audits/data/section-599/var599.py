# Audit 599: midpoint-orbit variants of rule B against exact QM, all E1 channels (hydrogen, atomic units).
import math
src=open('../section-560/tempo_e1_punkt_srodkowy.py').read().split('if __name__ == "__main__":')[0]
g={}; exec(src,g)
fc,rate_qm,channels,C,AU=g['fourier_components'],g['rate_qm'],g['channels'],g['C'],g['AU_RATE']
def nm_of(kind,n,n2):
    return {'A':0.5*(n+n2),'E':1/math.sqrt(0.5*(1/n**2+1/n2**2)),'H':1/(0.5*(1/n+1/n2)),'G':math.sqrt(n*n2)}[kind]
cache={}
def rate(kind,om,n,l,n2,l2):
    key=(kind,om,n,l,n2,l2)
    if key in cache: return cache[key]
    s=n-n2; nm=nm_of(kind,n,n2); Lm=max(l,l2)
    if Lm>=nm: cache[key]=0.0; return 0.0
    cp,cm=fc(nm,Lm,s,N=8192); c=cp if l2==l-1 else cm
    w=s/nm**3 if om=='o' else 0.5*(1/n2**2-1/n**2)
    r=2*max(l,l2)/(2*l+1)*2/3*c*c*w**3/C**3; cache[key]=r; return r
rows=[]
for kind in 'AEHG':
    for om in 'ox':
        lr=[math.log(rate(kind,om,n,l,n2,l2)/rate_qm(n,l,n2,l2)) for n in range(2,11) for l in range(n) for n2,l2 in channels(n,l) if rate(kind,om,n,l,n2,l2)>0]
        big=[rate(kind,om,n,1,1,0)/rate_qm(n,1,1,0) for n in range(3,9)]+[rate(kind,om,n,0,2,1)/rate_qm(n,0,2,1) for n in range(3,11)]
        lt=[]
        for n in range(2,9):
            for l in range(min(n,4)):
                if (n,l)==(2,0): continue
                q=sum(rate_qm(n,l,n2,l2) for n2,l2 in channels(n,l)); r=sum(rate(kind,om,n,l,n2,l2) for n2,l2 in channels(n,l))
                lt.append((q/r if r>0 else float('inf'),f"{n}{'spdf'[l]}"))
        worst=max(lt,key=lambda t:abs(t[0]-1))
        rows.append((kind+om,math.sqrt(sum(x*x for x in lr)/len(lr)),max(abs(x) for x in lr),min(big),max(big),worst))
        print(f"{kind}{om}: channels rms|ln| {rows[-1][1]:.3f} max|ln| {rows[-1][2]:.3f};  largest-Delta-n ratio {min(big):.3f}-{max(big):.3f};  lifetimes worst {worst[1]} tau/QM {worst[0]:.3f}",flush=True)
best=min(rows,key=lambda r:abs(r[5][0]-1)); print(f"\nbest one-place variant: {best[0]} (worst state {best[5][1]} {best[5][0]:.3f})")
