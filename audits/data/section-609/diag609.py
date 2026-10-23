# Audit 609: extend 605 -- Born weights from emission energy balance + rule-Q block weights + time averages of observables (QM diagonal ensemble vs classical TDVP).
import numpy as np, math
src=open('../section-601/qm601.py').read().split('print("B[T]')[0]; g={}; exec(src,g)
H2P,PT,J,Lv,S1,S2,Jsq=g['H2P'],g['PT'],g['J'],g['Lv'],g['S1'],g['S2'],g['Jsq']
S=[a+b for a,b in zip(S1,S2)]
def plus(X): return X[0]+1j*X[1]
def minus(X): return X[0]-1j*X[1]
herm=lambda A: A+A.conj().T
Proj=[np.outer(e,e) for e in np.eye(12)]
Ops=Proj+[herm(plus(S1)@minus(S2)), herm(plus(Lv)@minus(S)), herm(plus(Lv)@minus(S1))]
Mb=np.array([(1-i)+(0.5 if p==0 else -0.5)+(0.5 if q==0 else -0.5) for i in range(3) for p in range(2) for q in range(2)])
rng=np.random.default_rng(609)
def iso(n): v=rng.normal(size=(n,3)); return v/np.linalg.norm(v,axis=1,keepdims=True)
def lcoh(V):   # (N,3) -> (N,3) complex, l = 1 coherent along V
    c=V[:,2]; ph=np.arctan2(V[:,1],V[:,0]); s=np.sqrt(np.clip(1-c*c,0,None))
    return np.stack([np.exp(-1j*ph)*(1+c)/2, s/np.sqrt(2)+0j, np.exp(1j*ph)*(1-c)/2],1)
def spin(V):
    th=np.arccos(np.clip(V[:,2],-1,1)); ph=np.arctan2(V[:,1],V[:,0])
    return np.stack([np.cos(th/2)+0j, np.exp(1j*ph)*np.sin(th/2)],1)
def state(X):  # X (N,9): L^, n1, n2 (unit)
    a=lcoh(X[:,0:3]); b=spin(X[:,3:6]); c=spin(X[:,6:9])
    return np.einsum('ni,nj,nk->nijk',a,b,c).reshape(len(X),12)
def qsym(Op,X): psi=state(X); return np.einsum('ni,ij,nj->n',psi.conj(),Op,psi).real
def grad(H,X,eps=1e-6):
    G=np.zeros_like(X)
    for j in range(9):
        Xp=X.copy(); Xm=X.copy(); Xp[:,j]+=eps; Xm[:,j]-=eps
        for blk in (slice(0,3),slice(3,6),slice(6,9)):
            Xp[:,blk]/=np.linalg.norm(Xp[:,blk],axis=1,keepdims=True); Xm[:,blk]/=np.linalg.norm(Xm[:,blk],axis=1,keepdims=True)
        G[:,j]=(qsym(H,Xp)-qsym(H,Xm))/(2*eps)
    return G
def rhs(H,X):
    G=grad(H,X); out=np.zeros_like(X)
    for blk,jl in ((slice(0,3),1.0),(slice(3,6),0.5),(slice(6,9),0.5)):
        out[:,blk]=2*np.pi*np.cross(G[:,blk],X[:,blk])/jl
    return out
def solve(rows,rhsv,blocks):
    A=np.array(rows); y=np.array(rhsv); p=np.linalg.lstsq(A,y,rcond=None)[0]; p=np.clip(p,0,None)
    for idx,PM in blocks:
        s=p[idx].sum(); p[idx]=p[idx]*PM/s if s>0 else PM/len(idx)
    return p
print("B[T]  T1 max|p-Born|  T2: TV mean/max   3P p-Ps QM / TDVP   3P o-Ps QM / TDVP   (H_Q drift, J_z drift)")
for B in (0.003,0.01,0.03,0.1,0.3,1.0):
    H=H2P(B).real/1e9; w,V=np.linalg.eigh(H)          # GHz
    Mk=np.array([int(round(np.sum(np.abs(V[:,k])**2*Mb))) for k in range(12)])
    Odiag=[np.array([(V[:,k].conj()@O@V[:,k]).real for k in range(12)]) for O in Ops]
    PTk=np.array([(V[:,k].conj()@PT@V[:,k]).real for k in range(12)])
    X=np.concatenate([iso(60),iso(60),iso(60)],1); psi=state(X); born=np.abs(psi@V.conj())**2
    blocks=[]; rows0=[]
    for M in range(-2,3):
        idx=np.where(Mk==M)[0]; row=np.zeros(12); row[idx]=1; rows0.append(row); blocks.append(idx)
    # QM diagonal-ensemble averages (exact) and classical TDVP time averages
    H0=qsym(H,X); Jz0=X[:,2]*1.0+0.5*X[:,5]+0.5*X[:,8]
    dt=0.001 if B>=1.0 else 0.002; nT=int(10/dt); acc=[np.zeros(len(X)) for _ in Ops]; Y=X.copy()
    for _ in range(nT):
        k1=rhs(H,Y); k2=rhs(H,Y+0.5*dt*k1); k3=rhs(H,Y+0.5*dt*k2); k4=rhs(H,Y+dt*k3); Y=Y+dt/6*(k1+2*k2+2*k3+k4)
        for blk in (slice(0,3),slice(3,6),slice(6,9)): Y[:,blk]/=np.linalg.norm(Y[:,blk],axis=1,keepdims=True)
        for i,O in enumerate(Ops): acc[i]+=qsym(O,Y)
    tavg=[a/nT for a in acc]
    dH=np.max(np.abs(qsym(H,Y)-H0)); dJ=np.max(np.abs(Y[:,2]+0.5*Y[:,5]+0.5*Y[:,8]-Jz0))
    e1=0; tv=[]; pr=np.zeros(12); pc=np.zeros(12)
    for n in range(len(X)):
        PMs=[born[n][idx].sum() for idx in blocks]
        rows=list(rows0)+[w]+[Od for Od in Odiag]
        exact=[*PMs, np.sum(born[n]*w)]+[np.sum(born[n]*Od) for Od in Odiag]
        clas=[*PMs, H0[n]]+[t[n] for t in tavg]
        pq=solve(rows,exact,list(zip(blocks,PMs))); pcl=solve(rows,clas,list(zip(blocks,PMs)))
        e1=max(e1,np.max(np.abs(pq-born[n]))); tv.append(0.5*np.sum(np.abs(pcl-born[n])))
    # ensemble 3P shares from per-state singlet character (the configurations' weights are not needed for the ensemble: use T1-style per-state)
    sk=1-PTk; qp=np.sum(sk*(1-sk))/3; qo=np.sum((1-sk)**2)/9
    # TDVP-based ensemble: average reconstructed weights over configurations as estimate of p-Ps / o-Ps ensembles is not defined -> report TV only, plus
    # 3P share of the configuration ensemble (isotropic product states) Born vs reconstruction
    Pb=np.mean(born@PTk); Pc=np.mean([0]) 
    recon=[]
    for n in range(len(X)):
        PMs=[born[n][idx].sum() for idx in blocks]
        rows=list(rows0)+[w]+[Od for Od in Odiag]; clas=[*PMs, H0[n]]+[t[n] for t in tavg]
        recon.append(solve(rows,clas,list(zip(blocks,PMs)))@PTk)
    print(f"{B:5.3f}  {e1:.1e}        {np.mean(tv):.3f}/{np.max(tv):.3f}    3P (product-state ensemble) Born {Pb:.3f} / TDVP {np.mean(recon):.3f}    ({dH:.1e}, {dJ:.1e})",flush=True)
