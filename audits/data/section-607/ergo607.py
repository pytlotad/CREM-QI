# Audit 607: eigenstate ergodicity -- singlet character of 2P Ps field eigenstates from classical energy-shell slices (Weyl cells) vs QM.
import numpy as np
src=open('../section-601/qm601.py').read().split('print("B[T]')[0]; g={}; exec(src,g)
H2P,PT,J=g['H2P'],g['PT'],g['J']; I12=np.eye(12); PS=I12-PT
A=2.7415; Bt=0.45725; mu=9.284765e-24; h=6.62607015e-34
Mb=np.array([(1-i)+(0.5 if p==0 else -0.5)+(0.5 if q==0 else -0.5) for i in range(3) for p in range(2) for q in range(2)])
rng=np.random.default_rng(607); N=200000
def iso(n): v=rng.normal(size=(n,3)); return v/np.linalg.norm(v,axis=1,keepdims=True)
Lh,n1,n2=iso(N),iso(N),iso(N)
def Hcl(Llen,slen,a):
    L=Llen*Lh; s1=slen*n1; s2=slen*n2
    p1=np.sum(s1*Lh,1); p2=np.sum(s2*Lh,1); d12=np.sum(s1*s2,1)
    return A*np.sum(L*(s1+s2),1)+4*Bt*(1.5*(d12-p1*p2)-d12)+2*a*(s1[:,2]-s2[:,2]), L[:,2]+s1[:,2]+s2[:,2]
print("B[T]   3P p-Ps: QM   (i)    (ii)   |  3P o-Ps: QM   (i)    (ii)   | s_k QM vs (i), (ii) in block M=0 (energy order)")
for B in (0.003,0.01,0.03,0.1,0.3,1.0):
    a=mu*B/h/1e9
    w,V=np.linalg.eigh(H2P(B)); sk=np.array([(V[:,k].conj()@PS@V[:,k]).real for k in range(12)])
    Mk=np.array([int(round(np.sum(np.abs(V[:,k])**2*Mb))) for k in range(12)])
    qp=np.sum(sk*(1-sk))/3; qo=np.sum((1-sk)**2)/9
    # T1 check vs dephasing
    rp=PS/3; deph=sum(np.outer(V[:,k],V[:,k].conj())*(V[:,k].conj()@rp@V[:,k]) for k in range(12)); assert abs(np.trace(PT@deph).real-qp)<1e-12
    res={}
    for var,(Llen,slen) in (("i",(1.5,1.0)),("ii",(1.0,0.5))):
        E,Jz=Hcl(Llen,slen,a); Mc=np.clip(np.round(Jz),-2,2)
        if var=="i": sing=(np.linalg.norm(slen*n1+slen*n2,axis=1)<1.0).astype(float)
        else: sing=(1-np.sum(n1*n2,1))/4
        skc=np.zeros(12)
        for M in range(-2,3):
            ks=[k for k in range(12) if Mk[k]==M]; ks=sorted(ks,key=lambda k:w[k])
            sel=np.where(Mc==M)[0]; order=sel[np.argsort(E[sel])]; chunks=np.array_split(order,len(ks))
            for k,ch in zip(ks,chunks): skc[k]=np.mean(sing[ch])
        res[var]=(np.sum(skc*(1-skc))/3,np.sum((1-skc)**2)/9,skc)
    blk=sorted([k for k in range(12) if Mk[k]==0],key=lambda k:w[k])
    print(f"{B:5.3f}   {qp:.3f}  {res['i'][0]:.3f}  {res['ii'][0]:.3f}  |   {qo:.3f}  {res['i'][1]:.3f}  {res['ii'][1]:.3f}  | "
          +" ".join(f"{sk[k]:.2f}/{res['i'][2][k]:.2f}/{res['ii'][2][k]:.2f}" for k in blk))
