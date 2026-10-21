# Audit 546: Born projection of the coherent state of the classical configuration onto the field eigenstates of the model's 2P Hamiltonian (prereg546.txt).
import runpy, io, contextlib, numpy as np
buf=io.StringIO()
with contextlib.redirect_stdout(buf): ns=runpy.run_path('../section-542/qm542.py')
EJ,ES=ns['sets']['model']; PJ=ns['PJ']; PS=ns['PS']; PT=ns['PT']; h=ns['h']; g=ns['g']; muB=ns['muB']; Zop=ns['Zop']
H0=sum(e*1e9*h*P for e,P in zip(EJ,PJ))+ES*1e9*h*PS
rng=np.random.default_rng(546)
def coh_l1(u):   # |l=1, m=1 along u>: components m = +1, 0, -1
    th=np.arccos(np.clip(u[2],-1,1)); ph=np.arctan2(u[1],u[0])
    return np.array([np.exp(-1j*ph)*(1+np.cos(th))/2, np.sin(th)/np.sqrt(2), np.exp(1j*ph)*(1-np.cos(th))/2])
def coh_half(u):  # spin-1/2 up along u: components up, down
    th=np.arccos(np.clip(u[2],-1,1)); ph=np.arctan2(u[1],u[0])
    return np.array([np.exp(-1j*ph/2)*np.cos(th/2), np.exp(1j*ph/2)*np.sin(th/2)])
singlet=np.array([0,1,-1,0])/np.sqrt(2)
def iso(n): v=rng.normal(size=(n,3)); return v/np.linalg.norm(v,axis=1)[:,None]
def basinJ(l,s):  # 544 basin on J^ = L + hbar s^ (|L| = 1.5): rotate both
    J=1.5*l+s; Jh=J/np.linalg.norm(J); m=int(np.rint(1.5*Jh[2]))
    t=np.array([Jh[0],Jh[1],0.]) if m==0 else np.array([0,0,float(m)]); t/=np.linalg.norm(t)
    ax=np.cross(Jh,t); sa=np.linalg.norm(ax); ca=Jh@t
    if sa<1e-12: return l,s
    k=ax/sa; R=lambda v: v*ca+np.cross(k,v)*sa+k*(k@v)*(1-ca); return R(l),R(s)
N=20000; Ls=iso(N); Ss=iso(N)
for B in (0.01,0.03,0.1,0.3,1.0,3.0):
    H=H0+(g/2)*muB*B*Zop; w,V=np.linalg.eigh(H); sk=np.real(np.einsum('ik,ij,jk->k',V.conj(),PS,V))
    def PSof(psi): p=np.abs(V.conj().T@psi)**2; return p@sk, p.sum()
    trip=[]; sing=[]; bas=[]; perJ={0:[],1:[],2:[]}
    for i in range(N):
        l,s=Ls[i],Ss[i]
        trip.append(PSof(np.kron(coh_l1(l),np.kron(coh_half(s),coh_half(s))))[0])
        if i<4000:
            sing.append(PSof(np.kron(coh_l1(l),singlet))[0])
            lb,sb=basinJ(l,s); bas.append(PSof(np.kron(coh_l1(lb),np.kron(coh_half(sb),coh_half(sb))))[0])
            for J in (0,1,2):   # fixed |J_e| = J + 1/2: s^ at cos(S,L) from the triangle, azimuth random about l
                cb=((J+.5)**2-2.25-1)/3; e1=np.cross(l,[0.3,0.5,0.81]); e1/=np.linalg.norm(e1); e2=np.cross(l,e1); ph=rng.uniform(0,2*np.pi)
                sJ=cb*l+np.sqrt(max(0,1-cb*cb))*(np.cos(ph)*e1+np.sin(ph)*e2)
                perJ[J].append(PSof(np.kron(coh_l1(l),np.kron(coh_half(sJ),coh_half(sJ))))[0])
    r=[ns['PSdecay'](EJ,ES,B,j) for j in (0,1,2)]; qm=sum((2*j+1)*r[j] for j in (0,1,2))/9
    trip=np.array(trip); tot=trip.mean(); se=trip.std()/np.sqrt(N); pj=[np.mean(perJ[J]) for J in (0,1,2)]; mixJ=sum((2*J+1)*pj[J] for J in (0,1,2))/9
    print(f"B {B:5.2f} T: QM 3P {qm:.5e} | isotropic triplet {tot:.5e}+-{se:.1e} ({tot/qm:.4f}) | with 544 basin {np.mean(bas):.5e} ({np.mean(bas)/qm:.3f}) | "
          f"fixed |J_e| J0/J1/J2 {pj[0]/r[0]:.2f}/{pj[1]/r[1]:.2f}/{pj[2]/r[2]:.2f} x QM_J, (2J+1) mix {mixJ/qm:.3f} | singlet-born P_S {np.mean(sing):.4f} | "
          f"per-pair P_S spread (triplet) 10/50/90 % {np.percentile(trip,10):.4f}/{np.percentile(trip,50):.4f}/{np.percentile(trip,90):.4f}",flush=True)
