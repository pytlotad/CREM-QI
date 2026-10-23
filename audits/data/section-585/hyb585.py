# Audit 585: Delta l = +1 channels from l >= 1 added to the engine-style co-rotating count (432 E_all), offline, hydrogen.
import numpy as np
src=open('../section-560/tempo_e1_punkt_srodkowy.py').read().split('if __name__ == "__main__":')[0]
g={}; exec(src,g)
rate_qm, rate_rule, channels, C, AU = g['rate_qm'], g['rate_rule'], g['channels'], g['C'], g['AU_RATE']
def e_all(n,l,N=1<<20):
    a=n*n; L=l+0.5; e=np.sqrt(max(0.0,1-(L/n)**2)); b=a*np.sqrt(1-e*e)
    M=(np.arange(N))*2*np.pi/N; u=np.full_like(M,np.pi)
    for _ in range(200): u-= (u-e*np.sin(u)-M)/(1-e*np.cos(u))
    assert np.max(np.abs(u-e*np.sin(u)-M))<1e-12
    z=a*(np.cos(u)-e)+1j*b*np.sin(u)
    c=np.fft.fft(z)/N          # c[k] = <z e^{-ikM}> : co-rotating component of harmonic k
    k=np.arange(1,N//2); Om=1.0/n**3
    return np.sum(2/3*np.abs(c[k])**2*(k*Om)**3)/C**3
ns=lambda A: 1e9/(A*AU)
print("check E_all: 2p %.3f 3p %.3f 3d %.3f ns (432: 1.674 5.582 15.729)"%(ns(e_all(2,1)),ns(e_all(3,1)),ns(e_all(3,2))))
names="spdfghijk"
cE=cH=0; tot=0; maxH=0; brok=0; brtot=0
print("\nstate   tau_QM[ns]   E/QM    H/QM    F/QM   up-share QM   H    H/QM")
for n in range(3,11):
    for l in range(1,n-1):
        qd=sum(rate_qm(n,l,n2,l2) for n2,l2 in channels(n,l) if l2==l-1)
        qu=sum(rate_qm(n,l,n2,l2) for n2,l2 in channels(n,l) if l2==l+1)
        ru=sum(rate_rule(n,l,n2,l2,'RBw') for n2,l2 in channels(n,l) if l2==l+1)
        rf=sum(rate_rule(n,l,n2,l2,'RBw') for n2,l2 in channels(n,l))
        AE=e_all(n,l); q=qd+qu
        tE,tH,tF=q/AE,q/(AE+ru),q/rf
        su,sh=qu/q,ru/(AE+ru)
        tot+=1; cH+= abs(np.log(tH))<=abs(np.log(tE)); maxH=max(maxH,abs(tH-1))
        flag=''
        if su>=0.01:
            brtot+=1; ok=0.8<=sh/su<=1.25; brok+=ok; flag='' if ok else ' <-'
        print(f"{n}{names[l]}  {ns(q):10.3f}  {tE:6.3f}  {tH:6.3f}  {tF:6.3f}   {su:8.4f}  {sh:8.4f}  {sh/su:5.3f}{flag}")
print(f"\nH585a: H not worse than E in {cH}/{tot} ({100*cH/tot:.0f} %), max |tau_H/QM - 1| = {maxH:.3f}")
print(f"H585b: up-share H/QM in 0.8-1.25 in {brok}/{brtot} states with QM share >= 1 %")
