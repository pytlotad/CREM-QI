# Audit 602 (check after the run): conservation of H and J_z by the integrator of lyap602.py, 20 trajectories, 100 ns.
import numpy as np
src=open('lyap602.py').read().split('for slen,tag in')[0]; exec(src)
def H(Y,a):
    L,s1,s2=Y[:,0:3],Y[:,3:6],Y[:,6:9]; Lh=L/np.linalg.norm(L,axis=1,keepdims=True)
    p1=np.sum(s1*Lh,1); p2=np.sum(s2*Lh,1); d12=np.sum(s1*s2,1)
    return A*np.sum(L*(s1+s2),1)+4*Bt*(1.5*(d12-p1*p2)-d12)+2*a*(s1[:,2]-s2[:,2])
for slen in (1.0,0.5):
    for B in (0.0,0.1,1.0):
        rng=np.random.default_rng(1); a=mu*B/h/1e9; n=20
        Y=np.concatenate([1.5*iso(rng,n),slen*iso(rng,n),slen*iso(rng,n)],1); H0=H(Y,a); Jz0=Y[:,2]+Y[:,5]+Y[:,8]
        for _ in range(int(100/2e-3)): Y=rk4(Y,a,2e-3)
        print(f"|s|={slen} B={B}: max |dH| {np.max(np.abs(H(Y,a)-H0)):.1e} GHz (|H| ~ {np.max(np.abs(H0)):.1f}), max |dJz| {np.max(np.abs(Y[:,2]+Y[:,5]+Y[:,8]-Jz0)):.1e}",flush=True)
