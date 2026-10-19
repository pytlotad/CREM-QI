# Engine-like: isotropic L -> M uniform on [-(l+1/2), l+1/2]; allowed |M - dm| <= l - 1/2; weights: Lorentz at c = M/(l+1/2) (import 20 now)
# vs geometric mean at L_mid = l (clamped >= 0); basin m = round(M); per-basin mean of normalized weights vs CG^2.
import numpy as np
def cg(l,m,q):
    v={1:(l+m)*(l+m-1),0:2*(l*l-m*m),-1:(l-m)*(l-m-1)}[q]; return v/(2*l*(2*l-1))
for l in (2,3,4):
    M=np.linspace(-(l+0.5),l+0.5,400001)[1:-1]; Lp=l-0.5
    out={}
    for name in ('lorentz','geom'):
        W=[]
        for q in (1,0,-1):
            ok=np.abs(M-q)<=Lp+1e-12
            if name=='lorentz': c=M/(l+0.5); w=(1-c*c)/2 if q==0 else (1+q*c)**2/4
            else:
                ci=M/l; cf=(M-q)/l; w=(1-ci*ci)/2 if q==0 else (1+q*ci)*(1+q*cf)/4
            W.append(np.where(ok,np.maximum(w,0),0))
        W=np.array(W); S=W.sum(0); good=S>0; W=W/np.where(good,S,1)
        basin=np.rint(M).astype(int); worst=0; rows=[]
        for m in range(-l,l+1):
            sel=(basin==m)&good
            mean=W[:,sel].mean(1); ref=np.array([cg(l,m,q) if abs(m-q)<=l-1 else 0 for q in (1,0,-1)])
            worst=max(worst,np.abs(mean-ref).max()); rows.append((m,mean.round(4),ref.round(4)))
        tot=W[:,good].mean(1)
        out[name]=(worst,rows,tot,(~good).mean())
    for name,(worst,rows,tot,nog) in out.items():
        print(f"l {l} {name:8s}: max |basin mean - CG^2| {worst:.4f}; total sigma+/pi/sigma- {tot.round(4)} (QM 1/3 each); no allowed {nog:.4f}")
        if l==2: [print("     m",m,"model",a,"QM",b) for m,a,b in rows]
