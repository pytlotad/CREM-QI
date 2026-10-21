# Audit 545: eigenstates of the model's 2P Hamiltonian in a field, labelled by adiabatic continuation from B = 0 (energy order within each M block),
# triplet weights t_k, Born weights p_k = t_k/9; check sum p_k s_k against the QM sudden-start average of 542.  Writes eig545_<B>.txt for m545.cpp.
import runpy, io, contextlib, numpy as np
buf=io.StringIO()
with contextlib.redirect_stdout(buf): ns=runpy.run_path('../section-542/qm542.py')
EJ,ES=ns['sets']['model']; PJ=ns['PJ']; PS=ns['PS']; PT=ns['PT']; Jz=ns['J'][2].real; h=ns['h']; g=ns['g']; muB=ns['muB']; Zop=ns['Zop']
H0=sum(e*1e9*h*P for e,P in zip(EJ,PJ))+ES*1e9*h*PS
labels0=[]  # B = 0 labels per M block, energy order
for M in (-2,-1,0,1,2):
    lab=[('3P',J,EJ[J]) for J in (0,1,2) if abs(M)<=J]+([('1P',1,ES)] if abs(M)<=1 else [])
    labels0.append((M,sorted(lab,key=lambda x:x[2])))
for B in (0.01,0.03,0.1,0.3,1.0,3.0):
    H=H0+(g/2)*muB*B*Zop; out=[]; tot=0
    w,V=np.linalg.eigh(H+1e-9*h*Jz)   # tiny Jz term to make eigenvectors M-pure
    Mk=[np.real(V[:,i].conj()@Jz@V[:,i]) for i in range(12)]
    for M,lab in labels0:
        idx=sorted([i for i in range(12) if abs(Mk[i]-M)<1e-6],key=lambda i:w[i])
        assert len(idx)==len(lab)
        for i,(kind,J,_) in zip(idx,lab):
            t=np.real(V[:,i].conj()@PT@V[:,i]); out.append((kind,J,M,t,w[i]/h*1e-9)); tot+=t/9*(1-t)
    with open(f'eig545_{B:g}.txt','w') as f:
        for kind,J,M,t,E in out: f.write(f"{0 if kind=='3P' else 1} {J} {M} {t:.10f} {E:.6f}\n")
    r=[ns['PSdecay'](EJ,ES,B,j) for j in (0,1,2)]; qm=sum((2*j+1)*r[j] for j in (0,1,2))/9
    print(f"B {B:5.2f} T: sum (t_k/9) s_k = {tot:.5e}, QM sudden-start 3P {qm:.5e} (ratio {tot/qm:.4f}); "+' '.join(f"{k}{J}{M:+d}:t={t:.3f}" for k,J,M,t,E in out))
