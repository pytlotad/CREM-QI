# Audit 606: dephasing (projection) basis -- eigenbasis (Born) vs product (single-moment projections) vs coupled |S J M> (J basins).
import numpy as np
src=open('../section-601/qm601.py').read().split('print("B[T]')[0]; g={}; exec(src,g)
H2P,H1S,PT,Jsq,Ssq,J,trL,I3,sing,np_=g['H2P'],g['H1S'],g['PT'],g['Jsq'],g['Ssq'],g['J'],g['trL'],g['I3'],g['sing'],np
rho_p=np.kron(I3/3,np.outer(sing,sing)); PTspin=np.eye(4)-np.outer(sing,sing); rho_o=np.kron(I3/3,PTspin/3)
_,Vc=np.linalg.eigh(Jsq+0.0137*J[2]+10*Ssq)          # coupled basis |S J M>, nondegenerate labels
def dephase(rho,V): return sum(np.outer(V[:,k],V[:,k].conj())*(V[:,k].conj()@rho@V[:,k]) for k in range(V.shape[1]))
def tri1S(rho,B):
    r=trL(rho); w1,V1=np.linalg.eigh(H1S(B)); s=0
    for i in range(4):
        v=V1[:,i]
        if 0.5*abs(v[1]-v[2])**2>=0.5: s+=(v.conj()@r@v).real     # the mostly-singlet state S'
    return 1-s
print("B[T]   start   3P: EIG   PROD    CPL  |   1S triplet: EIG   PROD    CPL   | rule Q suffices?")
for B in (0.001,0.003,0.01,0.03,0.1,0.3,1.0,3.0,10.0,30.0):
    _,Ve=np.linalg.eigh(H2P(B)); Vp=np.eye(12)
    ok=True; rows=[]
    for name,r0 in (("p-Ps",rho_p),("o-Ps",rho_o)):
        vals={}
        for lab,V in (("EIG",Ve),("PROD",Vp),("CPL",Vc)):
            r=dephase(r0,V); vals[lab]=(np.trace(PT@r).real, tri1S(r,B))
        d=[min(abs(vals["PROD"][i]-vals["EIG"][i]),abs(vals["CPL"][i]-vals["EIG"][i])) for i in (0,1)]
        ok&= max(d)<=0.02
        rows.append(f"{B:6.3f}  {name}   {vals['EIG'][0]:.3f}  {vals['PROD'][0]:.3f}  {vals['CPL'][0]:.3f}  |   {vals['EIG'][1]:.3f}  {vals['PROD'][1]:.3f}  {vals['CPL'][1]:.3f}")
    print(rows[0]+f"   | {'TAK' if ok else 'NIE'}"); print(rows[1])
