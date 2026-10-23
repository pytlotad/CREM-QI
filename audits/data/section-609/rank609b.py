import numpy as np
src=open('diag609.py').read().split('print("B[T]  T1')[0]; g={}; exec(src,g)
H2P,Ops,Mb=g['H2P'],g['Ops'],g['Mb']
for B in (0.003,0.01,0.03,0.1,0.3,1.0):
    H=H2P(B).real/1e9; w,V=np.linalg.eigh(H); Mk=np.array([int(round(np.sum(np.abs(V[:,k])**2*Mb))) for k in range(12)])
    rows=[]
    for M in range(-2,3):
        r=np.zeros(12); r[Mk==M]=1; rows.append(r)
    rows.append(w); rows+=[np.array([(V[:,k].conj()@O@V[:,k]).real for k in range(12)]) for O in Ops]
    print(f"B={B}: rank {np.linalg.matrix_rank(np.array(rows),tol=1e-9)} of 12 ({len(rows)} equations)")
