# Audit 543: linear theory of the classical relative spin mode of Ps 2^3P_J in a weak field (prereg543.txt).
# P_S/QM = sum_k <W_k>/D_k^2 / (<CG^2>/Delta_J^2); W_0 = sin^2 b cos^2 a, W_+- = sin^2 a (1 +- cos b)^2/4; D_k = H - kappa w_c cos b + k w_c.
import numpy as np
E=np.array([-7.312222,-1.826996,2.558642]); X=np.array([0.5,1.5,2.5])
d01=(E[1]-E[0]); d12=(E[2]-E[1]); cu=(d12-d01)/2
slope=lambda x: d01+cu*((x-0.5)+(x-1.5))
rng=np.random.default_rng(543); N=400000
def ensemble(J,basin):
    cb=(1+(J+.5)**2-2.25)/(2*(J+.5)); sb=np.sqrt(max(0,1-cb*cb))   # cos(s^, J^)
    # random orientation: J^ direction u, s^ = cb u + sb (cos psi e1 + sin psi e2)
    u=rng.normal(size=(N,3)); u/=np.linalg.norm(u,axis=1)[:,None]
    t=np.cross(u,[0.3,0.5,0.81]); t/=np.linalg.norm(t,axis=1)[:,None]; t2=np.cross(u,t); psi=rng.uniform(0,2*np.pi,N)
    s=cb*u+sb*(np.cos(psi)[:,None]*t+np.sin(psi)[:,None]*t2)
    def rot_to(v,target,others):   # rotate rigid set so that v -> target (shortest arc)
        ax=np.cross(v,target); sa=np.linalg.norm(ax,axis=1); ca=np.sum(v*target,axis=1); k=np.where(sa[:,None]>1e-12,ax/np.maximum(sa,1e-300)[:,None],0)
        return [w*ca[:,None]+np.cross(k,w)*sa[:,None]+k*np.sum(k*w,axis=1)[:,None]*(1-ca)[:,None] for w in others]
    if basin in ('s','J'):
        ax=s if basin=='s' else u
        m=np.rint(1.5*ax[:,2]); tgt=np.where((m==0)[:,None],ax*[1,1,0],np.stack([0*m,0*m,np.sign(ax[:,2])],axis=1))
        tgt/=np.linalg.norm(tgt,axis=1)[:,None]; u,s=rot_to(ax,tgt,[u,s])
    ca=u[:,2]; sa2=1-ca*ca
    return cb, np.array([np.mean(sb*sb*ca*ca), np.mean(sa2*(1+cb)**2/4), np.mean(sa2*(1-cb)**2/4)])
for kappa,name in ((0,'C/U2 kinematics'),(1,'U3')):
    for basin in ('iso','s','J'):
        r=[]
        for J in (0,1,2):
            cb,W=ensemble(J,basin); H=E[J]; wc=slope(J+.5); Om=H-kappa*wc*cb; D=np.array([Om,Om+wc,Om-wc])
            r.append((np.sum(W/D**2)/((1/3)/H**2), W, D))
        qm=[(2*J+1)/9*(1/3)/E[J]**2 for J in (0,1,2)]; tot=sum(qm[J]*r[J][0] for J in (0,1,2))/sum(qm)
        print(f"{name:16s} basin {basin:3s}: J0 {r[0][0]:.3f} J1 {r[1][0]:.3f} J2 {r[2][0]:.3f} 3P {tot:.3f} | J1 W(k=0,+1,-1) {np.round(r[1][1],3)} D {np.round(r[1][2],3)} GHz")
