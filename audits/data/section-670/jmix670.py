# Audit 670: field at which the 2P eigenstates' j mixing max_k min(<P_3/2>, <P_1/2>) reaches the 1 % threshold (661 Hamiltonian, antihydrogen numbers
# of qm662.py: E_3/2 - E_1/2 10969 MHz, A_1/2 59.2, A_3/2 11.8 MHz, positron and orbit Zeeman -muB B, antiproton 21.4 MHz/T).  usage: python3 jmix670.py
import numpy as np
r=np.sqrt(.5)
lx=np.array([[0,r,0],[r,0,r],[0,r,0]]); ly=np.array([[0,-1j*r,0],[1j*r,0,-1j*r],[0,1j*r,0]]); lz=np.diag([1,0,-1.])
sx=np.array([[0,.5],[.5,0]]); sy=np.array([[0,-.5j],[.5j,0]]); sz=np.diag([.5,-.5]); I2=np.eye(2); I3=np.eye(3)
k3=lambda a,b,c: np.kron(np.kron(a,b),c)
Lv=[k3(o,I2,I2) for o in (lx,ly,lz)]; Iv=[k3(I3,o,I2) for o in (sx,sy,sz)]; Sv=[k3(I3,I2,o) for o in (sx,sy,sz)]
J=[Lv[i]+Sv[i] for i in range(3)]; J2=sum(j@j for j in J); JI=sum(J[i]@Iv[i] for i in range(3)); Id=np.eye(12)
P1=(J2-3.75*Id)/(-3); P3=(J2-.75*Id)/3
def jmix(B):
    muB=13996.245*B; H=10969.0*P3+59.2*P1@JI+11.8*P3@JI+2*21.4*B*Iv[2]+2*(-1.00116*muB)*Sv[2]+(-muB)*Lv[2]
    w,V=np.linalg.eigh(H); p=np.real(np.einsum('ik,ij,jk->k',V.conj(),P3,V)); return np.max(np.minimum(p,1-p))
for B in (1e-6,1e-3,0.01,0.03,0.05,0.1,0.3,1.0329): print(f'B {B:9.4g} T: jmix {jmix(B):.3e}')
lo,hi=1e-4,1.0
for _ in range(60):
    m=np.sqrt(lo*hi)
    if jmix(m)<0.01: lo=m
    else: hi=m
print(f'threshold 1 %: B = {lo:.4f} T')
