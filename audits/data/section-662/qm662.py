# Audit 662: 1S hyperfine populations after 2P -> 1S from |L^> (x) |S> (spin singlet), model's 2P Hamiltonian (audit 661 numbers, H at 1e-6 T ~ 0 field):
# (a) the rule of 661/662 -- Born projection on 2P eigenstates, 1S from Tr_L |k><k| (full dephasing);  (b) QM reference: coherent decay with Gamma
# (tau 2P = 1.596 ns measured), the k != k' terms weighted Gamma/(Gamma + i omega_kk').  usage: python3 qm662.py [B tesla]
import numpy as np,sys
B=float(sys.argv[1]) if len(sys.argv)>1 else 0.0
MHz=1.0  # energies in MHz
E12,E32=0.0,10969.0; A12,A32=59.2,11.8; muB=13996.245*B; aL=-muB; ae=-1.00116*muB; ap=0.0   # antihydrogen signs as 661 (aL = a_e+ ~ -muB B)
ap=21.4*B
K=1420.4/4; tau=1.596e-9; G=1/(2*np.pi*tau)/1e6   # Gamma/2pi in MHz
r=np.sqrt(.5)
lx=np.array([[0,r,0],[r,0,r],[0,r,0]]); ly=np.array([[0,-1j*r,0],[1j*r,0,-1j*r],[0,1j*r,0]]); lz=np.diag([1,0,-1.])
sx=np.array([[0,.5],[.5,0]]); sy=np.array([[0,-.5j],[.5j,0]]); sz=np.diag([.5,-.5]); I2=np.eye(2); I3=np.eye(3)
k3=lambda a,b,c: np.kron(np.kron(a,b),c)
# first = antiproton (heavy), second = positron (light), as PAIR=antiproton,positron
Lv=[k3(o,I2,I2) for o in (lx,ly,lz)]; Iv=[k3(I3,o,I2) for o in (sx,sy,sz)]; Sv=[k3(I3,I2,o) for o in (sx,sy,sz)]
J=[Lv[i]+Sv[i] for i in range(3)]; J2=sum(j@j for j in J); JI=sum(J[i]@Iv[i] for i in range(3)); Id=np.eye(12)
P1=(J2-3.75*Id)/(-3); P3=(J2-.75*Id)/3
H=E12*P1+E32*P3+A12*P1@JI+A32*P3@JI+2*ap*Iv[2]+2*ae*Sv[2]+aL*Lv[2]
w,V=np.linalg.eigh(H)
singlet=np.array([0,r,-r,0])  # |m_pbar m_e+>: up up, up down, down up, down down
H1=K*4*sum(np.kron(a,b) for a,b in ((sx,sx),(sy,sy),(sz,sz)))+2*ap*np.kron(sz,I2)+2*ae*np.kron(I2,sz)
w1,V1=np.linalg.eigh(H1)
def lstate(th,ph):
    c=np.cos(th); e=np.exp(-1j*ph); return np.array([e*(1+c)/2,np.sqrt(1-c*c)/np.sqrt(2),np.conj(e)*(1-c)/2])
def run(th,ph):
    psi=np.kron(lstate(th,ph),singlet); c=V.conj().T@psi
    rule=np.zeros((4,4),complex); coh=np.zeros((4,4),complex)
    for k in range(12):
        for kk in range(12):
            vk=V[:,k].reshape(3,4); vkk=V[:,kk].reshape(3,4)
            T=sum(np.outer(vk[m],vkk[m].conj()) for m in range(3))   # Tr_L |k><kk|
            term=c[k]*np.conj(c[kk])*T
            coh+=term*G/(G+1j*(w[k]-w[kk]))
            if k==kk: rule+=term
    pop=lambda rho: np.real(np.einsum('ai,ab,bi->i',V1.conj(),rho,V1))
    return pop(rule),pop(coh)
th=np.arccos(np.linspace(-0.95,0.95,20)); R=np.zeros(4); C=np.zeros(4)
for t in th:
    for p in (0.3,1.9,4.1):
        a,b=run(t,p); R+=a; C+=b
R/=R.sum(); C/=C.sum()
print(f'B = {B} T; 1S eigenstates (MHz rel.): '+' '.join(f'{x:.1f}' for x in w1))
print('  singlet weight of each 1S eigenstate: '+' '.join(f'{abs(V1[:,i]@singlet)**2:.3f}' for i in range(4)))
print('  (a) rule (eigen projection) populations: '+' '.join(f'{x:.4f}' for x in R))
print('  (b) QM coherent decay, Gamma/2pi = %.1f MHz: '%G+' '.join(f'{x:.4f}' for x in C))
