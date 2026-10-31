# Audit 667: projection inside the classical j manifold (rule R1) for |L^> (x) |S>, model's 2P Hamiltonian at zero field (661 numbers, as qm662.py).
# usage: python3 qm667.py <662 b0 dir>
import numpy as np,sys,re,glob,math
E12,E32=0.0,10969.0; A12,A32=59.2,11.8; K=1420.4/4      # MHz
r=np.sqrt(.5)
lx=np.array([[0,r,0],[r,0,r],[0,r,0]]); ly=np.array([[0,-1j*r,0],[1j*r,0,-1j*r],[0,1j*r,0]]); lz=np.diag([1,0,-1.])
sx=np.array([[0,.5],[.5,0]]); sy=np.array([[0,-.5j],[.5j,0]]); sz=np.diag([.5,-.5]); I2=np.eye(2); I3=np.eye(3)
k3=lambda a,b,c: np.kron(np.kron(a,b),c)
Lv=[k3(o,I2,I2) for o in (lx,ly,lz)]; Iv=[k3(I3,o,I2) for o in (sx,sy,sz)]; Sv=[k3(I3,I2,o) for o in (sx,sy,sz)]
J=[Lv[i]+Sv[i] for i in range(3)]; J2=sum(j@j for j in J); JI=sum(J[i]@Iv[i] for i in range(3)); Id=np.eye(12)
P1=(J2-3.75*Id)/(-3); P3=(J2-.75*Id)/3
H=E12*P1+E32*P3+A12*P1@JI+A32*P3@JI
w,V=np.linalg.eigh(H)
jlab=np.array([1.5 if np.real(V[:,k].conj()@P3@V[:,k])>0.5 else 0.5 for k in range(12)])
singlet=np.array([0,r,-r,0])
H1=K*4*sum(np.kron(a,b) for a,b in ((sx,sx),(sy,sy),(sz,sz))); w1,V1=np.linalg.eigh(H1)
def lstate(th,ph):
    c=np.cos(th); e=np.exp(-1j*ph); return np.array([e*(1+c)/2,np.sqrt(1-c*c)/np.sqrt(2),np.conj(e)*(1-c)/2])
pop=lambda rho: np.real(np.einsum('ai,ab,bi->i',V1.conj(),rho,V1))
def rule(psi,mask):
    c=V.conj().T@psi; rho=np.zeros((4,4),complex)
    for k in range(12):
        if not mask[k]: continue
        vk=V[:,k].reshape(3,4); rho+=abs(c[k])**2*sum(np.outer(vk[m],vk[m].conj()) for m in range(3))
    return pop(rho)
th=np.arccos(np.linspace(-0.95,0.95,20)); Pj={0.5:0.0,1.5:0.0}; F={0.5:np.zeros(4),1.5:np.zeros(4)}; Rall=np.zeros(4); N=0
for t in th:
    for p in (0.3,1.9,4.1):
        psi=np.kron(lstate(t,p),singlet); N+=1
        for j in (0.5,1.5):
            Pj[j]+=np.real(psi.conj()@(P1 if j==0.5 else P3)@psi); F[j]+=rule(psi,jlab==j)
        Rall+=rule(psi,np.ones(12,bool))
for j in Pj: Pj[j]/=N
f={j:F[j][0]/F[j].sum() for j in F}      # 1S state 0 = F = 0 (singlet weight 1 at zero field)
tot=sum(Pj[j]*f[j] for j in Pj)
print(f'STEP A: <P_1/2> {Pj[0.5]:.4f}  <P_3/2> {Pj[1.5]:.4f};  f_1/2 (F=0 share from j=1/2) {f[0.5]:.4f}  f_3/2 {f[1.5]:.4f}')
print(f'  sum_j <P_j> f_j = {tot:.4f}  (662 rule, all j: {Rall[0]/Rall.sum():.4f}; QM 5/9 = 0.5556): H667a {abs(tot-5/9)<=1e-3}')
for j in (0.5,1.5):
    ww=w[jlab==j]; print(f'  j={j}: eigenvalues {np.round(ww,2)} MHz, spread {ww.max()-ww.min():.2f} MHz = {(ww.max()-ww.min())*4.135667696e-3:.4f} ueV')
spread=max(w[jlab==j].max()-w[jlab==j].min() for j in (0.5,1.5))*4.135667696e-3
print(f'  H667d: R1 import at zero field <= spread {spread:.4f} ueV <= 1 ueV: {spread<=1}')
# STEP B: engine's classical j from 662 b0
u=1.602176634e-25; n32=0; n=0
for fn in sorted(glob.glob(sys.argv[1]+'/A1_*.err')):
    t=open(fn).read()
    d=re.search(r'PAIR eigen draw k=\d+ P=\S+ dm \S+ E (\S+) (\S+) .* w_k (\S+)',t); im=re.search(r'PAIR eigen photon .* import (\S+)',t)
    if not (d and im): continue
    E0,E1,wk=float(d.group(1)),float(d.group(2)),float(d.group(3)); imp=float(im.group(1))/u
    jdraw=0.5 if abs(wk-E0)<abs(wk-E1) else 1.5
    jcl=jdraw if abs(imp)<1.0 else (2.0-jdraw)
    n+=1; n32+=jcl==1.5
p=n32/n; s=math.sqrt((2/3)*(1/3)/n)
print(f'STEP B: engine classical j = 3/2 in {n32}/{n} = {p:.3f} (QM <P_3/2> 0.667 +- {2*s:.3f} at 2 sigma): H667b {abs(p-2/3)<=2*s}')
pred=(1-p)*f[0.5]+p*f[1.5]; sp=math.sqrt(p*(1-p)/n)*abs(f[1.5]-f[0.5])
print(f'  R1 with the engine frequencies: F=0 share {pred:.4f} +- {sp:.4f} (QM 0.5556): H667c {abs(pred-5/9)<=2*sp}')
