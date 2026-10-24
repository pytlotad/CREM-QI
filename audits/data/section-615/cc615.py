# Audit 615: the author's principle |c|^2 -> c_i c_f (amplitudes at both ends of the transition), lifetimes vs exact QM (hydrogen).
import math
src=open('../section-560/tempo_e1_punkt_srodkowy.py').read().split('if __name__ == "__main__":')[0]; g={}; exec(src,g)
fc,rate_qm,channels,C=g['fourier_components'],g['rate_qm'],g['channels'],g['C']
def amp(n,L,s,co):
    cp,cm=fc(n,L,s,N=8192); return cp if co else cm
def rate(v,n,l,n2,l2):
    s=n-n2; co=(l2==l-1); w=2*max(l,l2)/(2*l+1)
    ci=amp(n,l+0.5,s,co); cf=amp(n2,l2+0.5,s,co)
    om=0.5*(1/n2**2-1/n**2) if v=='Bohr' else s/(0.5*(n+n2))**3
    return w*2/3*ci*cf*om**3/C**3
for v in ('Bohr','mid'):
    out=[]
    for name,(n,l) in (("2p",(2,1)),("3s",(3,0)),("3p",(3,1)),("4p",(4,1)),("3d",(3,2)),("5s",(5,0)),("8p",(8,1)),("8d",(8,2))):
        q=sum(rate_qm(n,l,*c) for c in channels(n,l)); r=sum(rate(v,n,l,*c) for c in channels(n,l)); out.append(f"{name} {q/r:.3f}")
    print(f"omega {v}: tau/QM "+"  ".join(out))
