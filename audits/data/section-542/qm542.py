# Audit 542: QM reference, Ps n = 2, l = 1 in a field: 12 states |m_L>|m1 m2>, H = sum_J E_J P(3P_J) + E_S P(1P_1) + Z (sigma1z - sigma2z).
# Singlet fraction at the 2P decay, P_S = Gamma int e^{-Gamma t} Tr[P_S rho(t)] dt, unpolarized 3P_J start (rho0 = P_J/(2J+1)).
import numpy as np
h=6.62607015e-34; hbar=h/(2*np.pi); muB=9.2740100783e-24; g=2.0023193044; tau=3.348e-9; G=1/tau
sx=np.array([[0,1],[1,0]])/2; sy=np.array([[0,-1j],[1j,0]])/2; sz=np.diag([0.5,-0.5]); I2=np.eye(2)
s2=1/np.sqrt(2); Lx=np.array([[0,s2,0],[s2,0,s2],[0,s2,0]]); Ly=np.array([[0,-1j*s2,0],[1j*s2,0,-1j*s2],[0,1j*s2,0]]); Lz=np.diag([1.,0,-1]); I3=np.eye(3)
k=lambda a,b,c_: np.kron(np.kron(a,b),c_)
L=[k(Lx,I2,I2),k(Ly,I2,I2),k(Lz,I2,I2)]; S1=[k(I3,sx,I2),k(I3,sy,I2),k(I3,sz,I2)]; S2=[k(I3,I2,sx),k(I3,I2,sy),k(I3,I2,sz)]
S=[a+b for a,b in zip(S1,S2)]; J=[a+b for a,b in zip(L,S)]
Ssq=sum(x@x for x in S).real; Jsq=sum(x@x for x in J).real
PT=Ssq/2; PS=np.eye(12)-PT
def proj(Jv):  # projector onto 3P_J
    P=PT.copy()
    for Jo in (0,1,2):
        if Jo!=Jv: P=P@(Jsq-Jo*(Jo+1)*np.eye(12))/(Jv*(Jv+1)-Jo*(Jo+1))
    return P@PT
PJ=[proj(j) for j in (0,1,2)]
for j in (0,1,2): assert abs(np.trace(PJ[j]).real-(2*j+1))<1e-9
sets={'model':([-7.312222,-1.826996,2.558642],0.0),'measured':([-18.500+11.180,-13.012+11.180,-8.624+11.180],0.0)}
Zop=2*(S1[2]-S2[2])
def PSdecay(EJ,ES,B,j):
    H=sum(e*1e9*h*P for e,P in zip(EJ,PJ))+ES*1e9*h*PS+(g/2)*muB*B*Zop
    w,V=np.linalg.eigh(H); rho0=V.conj().T@PJ[j]@V/(2*j+1); Pe=V.conj().T@PS@V
    om=(w[:,None]-w[None,:])/hbar; rbar=rho0*G/(G+1j*om)
    return np.real(np.sum(rbar*Pe.T))
Bs=(0.0,0.01,0.03,0.1,0.3,1.0,3.0)
for name,(EJ,ES) in sets.items():
    print(f"QM levels {name}: E_J (J=0,1,2) {EJ} GHz rel. 1P_1, tau(2P) {tau*1e9} ns")
    for B in Bs:
        r=[PSdecay(EJ,ES,B,j) for j in (0,1,2)]; tot=sum((2*j+1)*r[j] for j in (0,1,2))/9
        Z=(g/2)*muB*B; wf=[(8/3)*(Z/(abs(EJ[j]-ES)*1e9*h))**2 for j in (0,1,2)]
        print(f"  B {B:5.2f} T: P_S J0 {r[0]:.5e} J1 {r[1]:.5e} J2 {r[2]:.5e} 3P {tot:.5e} | weak-field (8/3)(Z/D)^2 {wf[0]:.3e} {wf[1]:.3e} {wf[2]:.3e}")
