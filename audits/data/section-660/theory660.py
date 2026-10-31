# hydrogen-like 1S -> 2P lines at B (no Lamb, no 2P HFS), positron sign convention via CPT = hydrogen; energies relative to Dirac zero-field 1S(centroid)->2P1/2
import numpy as np
B=1.0329; muB=13.996245e3*B  # MHz
g=2.00231930; FS=10969.045; a=2*FS/3; A1=1420.4057
# 1S states at high field (electron m_s), with HFS diagonal +-A/4 approx for |up up> d; use exact Breit-Rabi for d and a
def S1(ms):  # d: ms=+1/2,mI=+1/2 pure: +A/4 + g muB B/2 ; a(ms=-1/2,mI=-1/2) pure: +A/4 - g muB B/2
    return A1/4 + ms*g*muB
lines=[]
for mj in (1.5,0.5,-0.5,-1.5):
    basis=[(ml,ms) for ml in (-1,0,1) for ms in (0.5,-0.5) if ml+ms==mj]
    H=np.zeros((len(basis),)*2)
    for i,(l1,s1) in enumerate(basis):
        for j,(l2,s2) in enumerate(basis):
            if i==j: H[i,i]=a*l1*s1 + muB*(l1+g*s1)
            elif l2==l1+1 and s2==s1-1: H[i,j]=a*0.5*np.sqrt(2-l1*(l1+1))
            elif l2==l1-1 and s2==s1+1: H[i,j]=a*0.5*np.sqrt(2-l1*(l1-1))
    w,v=np.linalg.eigh(H)
    for k in range(len(w)):
        for ms in (0.5,-0.5):
            for (ml,s),c in zip(basis,v[:,k]):
                if s==ms and c*c>1e-6:
                    lines.append((w[k]+a-S1(ms) , mj, ms, ml, c*c))  # +a: 2P1/2 zero-field at 0
for L in sorted(lines): print(f'{L[0]:10.1f} MHz  mj {L[1]:+.1f}  1S ms {L[2]:+.1f}  m_l {L[3]:+d} (dm photon {L[3]:+d}) weight {L[4]:.3f}')
