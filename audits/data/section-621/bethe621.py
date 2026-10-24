# Audit 621: Bethe logarithm of hydrogen s states from the classical orbit spectrum (V1, V2) and from rule-B differences (V3) vs exact ln k0.
import numpy as np, math
src=open('../section-560/tempo_e1_punkt_srodkowy.py').read().split('if __name__ == "__main__":')[0]; g={}; exec(src,g)
fc=g['fourier_components']
exact={1:2.984128556,2:2.811769893,3:2.767663612,4:2.749811840}
alpha=7.2973525693e-3; Ry=3.2898419603e15
def spectrum(n,L=0.5,N=1<<20):
    a=n*n; e=math.sqrt(max(0.0,1-(L/n)**2)); b=a*math.sqrt(1-e*e)
    M=np.arange(N)*2*np.pi/N; u=np.full_like(M,np.pi)
    for _ in range(200): u-=(u-e*np.sin(u)-M)/(1-e*np.cos(u))
    z=a*(np.cos(u)-e)+1j*b*np.sin(u); c=np.fft.fft(z)/N; k=np.arange(1,N//2)
    return k, np.abs(c[k])**2+np.abs(c[N-k])**2           # co + counter (s state: both modes)
res={}
for n in (1,2,3,4):
    k,c2=spectrum(n); E=2*k/n**3                          # Ry
    v1=np.sum(c2*k**4*np.log(E))/np.sum(c2*k**4)
    v2=np.sum(c2*k**3*np.log(E))/np.sum(c2*k**3)
    num=den=0.0
    for n2 in range(n+1,401):                             # absorption ns -> n'p
        s=n2-n; nm=0.5*(n+n2); cp,cm=fc(nm,1.0,s,N=1<<14); dE=1/n**2-1/n2**2; w=2*cp*cp*dE**3
        num+=w*math.log(dE); den+=w
    for n2 in range(2,n):                                 # emission ns -> n'p
        s=n-n2; nm=0.5*(n+n2); cp,cm=fc(nm,1.0,s,N=1<<14); dE=1/n2**2-1/n**2; w=2*cm*cm*dE**3
        num-=w*math.log(dE); den-=w
    v3=num/den
    res[n]=(v1,v2,v3)
    print(f"{n}s: exact {exact[n]:.4f}  V1 {v1:.4f} ({v1/exact[n]-1:+.1%})  V2 {v2:.4f} ({v2/exact[n]-1:+.1%})  V3 {v3:.4f} ({v3/exact[n]-1:+.1%})",flush=True)
pre=alpha**3/(3*math.pi)*Ry/1e6                           # MHz per unit of ln k0 at n = 2
print(f"\nLamb 2S: (alpha^3/3pi) Ry = {pre:.1f} MHz per unit ln k0; leading 2S with exact ln k0: {pre*(math.log(1/alpha**2)-exact[2]+19/30):.1f} MHz (measured 2S-2P 1057,845(9) MHz incl. 2P and higher orders)")
for i,v in enumerate(("V1","V2","V3")):
    d=res[2][i]-exact[2]; print(f"  {v}: ln k0(2s) error {d:+.3f} -> 2S shift error {-pre*d:+.1f} MHz ({-pre*d/1057.845:+.1%} of the measured interval)")
