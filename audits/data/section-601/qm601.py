# Audit 601: QM shares of the 1S field eigenstates after the 2P -> 1S photon, start p-Ps 2^1P with isotropic L (rho0 = I_L/3 x |S><S|).
import numpy as np, math
h=6.62607015e-34
EJ=np.array([-7.312222,-1.826996,2.558642])*1e9          # Hz, from p601.txt
dE=203.3811e9; mu=9.284765e-24
Gamma=1/3.348e-9                                          # 1/s, engine photon time 2p -> 1s (589)
r=math.sqrt(0.5)
lx=np.array([[0,r,0],[r,0,r],[0,r,0]]); ly=np.array([[0,-1j*r,0],[1j*r,0,-1j*r],[0,1j*r,0]]); lz=np.diag([1.,0,-1])
sx=np.array([[0,.5],[.5,0]]); sy=np.array([[0,-.5j],[.5j,0]]); sz=np.diag([.5,-.5]); I2=np.eye(2); I3=np.eye(3)
k3=lambda a,b,c: np.kron(np.kron(a,b),c)
S1=[k3(I3,s,I2) for s in (sx,sy,sz)]; S2=[k3(I3,I2,s) for s in (sx,sy,sz)]; Lv=[k3(l,I2,I2) for l in (lx,ly,lz)]
S=[a+b for a,b in zip(S1,S2)]; J=[a+b for a,b in zip(Lv,S)]
Ssq=sum(x@x for x in S); Jsq=sum(x@x for x in J); PT=Ssq/2; I12=np.eye(12)
def H2P(B):
    H=np.zeros((12,12),complex)
    for j in range(3):
        P=PT.copy()
        for K in range(3):
            if K!=j: P=P@(Jsq-K*(K+1)*I12)/(j*(j+1)-K*(K+1))
        H+=EJ[j]*P
    a=mu*B/h                                              # a1 = +mu B (electron), a2 = -mu B, per sigma_z
    return H+2*a*S1[2]-2*a*S2[2]
def H1S(B):
    K=dE/4; a=mu*B/h
    return np.array([[K+0,0,0,0],[0,-K+2*a,2*K,0],[0,2*K,-K-2*a,0],[0,0,0,K]],complex)  # a1 + a2 = 0, a1 - a2 = 2a
sing=np.zeros(4); sing[1]=r; sing[2]=-r
rho0=np.kron(I3/3,np.outer(sing,sing))
def trL(M): return sum(M[4*i:4*i+4,4*i:4*i+4] for i in range(3))
print("B[T]   secular: T+1 T-1 T0' S'  (triplet)        QM coherent: T+1 T-1 T0' S'  (triplet)   |diff triplet|")
for B in (0.01,0.03,0.1,0.3,1.0,3.0):
    w,V=np.linalg.eigh(H2P(B)); R=V.conj().T@rho0@V
    sec=sum(R[k,k].real*trL(np.outer(V[:,k],V[:,k].conj())) for k in range(12))
    coh=np.zeros((4,4),complex)
    for k in range(12):
        for kp in range(12):
            if abs(R[k,kp])<1e-15: continue
            fac=Gamma/(Gamma+1j*2*np.pi*(w[k]-w[kp]))
            coh+=R[k,kp]*fac*trL(np.outer(V[:,k],V[:,kp].conj()))
    w1,V1=np.linalg.eigh(H1S(B))
    def shares(rho):
        out={'T+1':0,'T-1':0,"T0'":0,"S'":0}
        for i in range(4):
            v=V1[:,i]; q=(v.conj()@rho@v).real
            if abs(v[0])**2>0.5: out['T+1']+=q
            elif abs(v[3])**2>0.5: out['T-1']+=q
            elif 0.5*abs(v[1]-v[2])**2>=0.5: out["S'"]+=q
            else: out["T0'"]+=q
        return out
    a=shares(sec); b=shares(coh); ta=1-a["S'"]; tb=1-b["S'"]
    print(f"{B:5.2f}  "+" ".join(f"{a[x]:.5f}" for x in a)+f"  ({ta:.5f})    "+" ".join(f"{b[x]:.5f}" for x in b)+f"  ({tb:.5f})   {abs(ta-tb):.1e}")
