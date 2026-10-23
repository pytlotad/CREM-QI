# Audit 609 (check after T1 failed): rank of the equation set, and of the alternative set of 12 product-basis projectors.
import numpy as np
src=open('diag609.py').read().split('print("B[T]  T1')[0]; g={}; exec(src,g)
H2P,Ops,Mb=g['H2P'],g['Ops'],g['Mb']
for B in (0.003,0.03,0.1,0.3,1.0):
    H=H2P(B).real/1e9; w,V=np.linalg.eigh(H)
    Mk=np.array([int(round(np.sum(np.abs(V[:,k])**2*Mb))) for k in range(12)])
    rows=[]
    for M in range(-2,3):
        r=np.zeros(12); r[Mk==M]=1; rows.append(r)
    rows.append(w); rows+=[np.array([(V[:,k].conj()@O@V[:,k]).real for k in range(12)]) for O in Ops]
    A=np.array(rows); P=np.abs(V)**2
    print(f"B={B}: 12-equation set rank {np.linalg.matrix_rank(A,tol=1e-9)}; 12 product projectors rank {np.linalg.matrix_rank(P,tol=1e-9)}, cond {np.linalg.cond(P):.1e}")
