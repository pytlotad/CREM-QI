# Audit 669: field rules for the 2P -> 1S projection against the coherent-decay QM from the same classical states (offline, from 668 probes).
# usage: python3 qm669.py <668 runs dir>
import numpy as np,re,glob,sys,math
u=1.602176634e-25; hP=6.62607015e-34; G=99.7                   # ueV in J; Gamma/2pi in MHz
r=np.sqrt(.5)
lx=np.array([[0,r,0],[r,0,r],[0,r,0]],complex); ly=np.array([[0,-1j*r,0],[1j*r,0,-1j*r],[0,1j*r,0]]); lz=np.diag([1,0,-1.]).astype(complex)
sx=np.array([[0,.5],[.5,0]],complex); sy=np.array([[0,-.5j],[.5j,0]]); sz=np.diag([.5,-.5]).astype(complex); I2=np.eye(2); I3=np.eye(3)
k3=lambda a,b,c: np.kron(np.kron(a,b),c)
Lv=[k3(o,I2,I2) for o in (lx,ly,lz)]; S1=[k3(I3,o,I2) for o in (sx,sy,sz)]; S2=[k3(I3,I2,o) for o in (sx,sy,sz)]
J=[Lv[i]+S2[i] for i in range(3)]; J2=sum(j@j for j in J); JI=sum(J[i]@S1[i] for i in range(3)); Id=np.eye(12)
P1=(J2-3.75*Id)/(-3); P3=(J2-.75*Id)/3
Mj=np.array([(1-j//4)+(0.5-(j//2)%2)+(0.5-j%2) for j in range(12)])
s1v=[np.kron(o,I2) for o in (sx,sy,sz)]; s2v=[np.kron(I2,o) for o in (sx,sy,sz)]; SS=sum(s1v[i]@s2v[i] for i in range(3))
singlet=np.array([0,r,-r,0])
def v3(t,key,n=3): m=re.search(key+r' '+' '.join([r'(\S+)']*n),t); return np.array([float(x) for x in m.groups()])
def coherentL(uv):
    ct=np.clip(uv[2],-1,1); ph=math.atan2(uv[1],uv[0]); e=np.exp(-1j*ph)
    return np.array([e*(1+ct)/2,math.sqrt(max(0,1-ct*ct))/math.sqrt(2),np.conj(e)*(1-ct)/2])
def coherentHalf(uv):
    th=math.acos(np.clip(uv[2],-1,1)); ph=math.atan2(uv[1],uv[0])
    return np.array([math.cos(th/2)*np.exp(-0.5j*ph),math.sin(th/2)*np.exp(0.5j*ph)])
def frame(bv):
    b=bv/np.linalg.norm(bv); t=np.array([1.,0,0]) if abs(b[0])<0.9 else np.array([0,1.,0])
    e1=t-b*(t@b); e1/=np.linalg.norm(e1); return e1,np.cross(b,e1),b
res={1:[],2:[]}; val=[]
for fn in sorted(glob.glob(sys.argv[1]+'/A*_*.err')):
    t=open(fn).read(); ph=int(re.search(r'/A(\d)_',fn).group(1))
    d=re.search(r'PAIR eigen draw k=\d+ P=\S+ dm (\S+) E (\S+) (\S+) A (\S+) (\S+) aL (\S+) a1 (\S+) a2 (\S+) K (\S+)',t)
    pe=re.search(r'PAIR eigen photon .* zeemanOrbit (\S+) import (\S+) gapShift (\S+)',t)
    lg=re.search(r'PAIR eigen ledger Hbar (\S+) Hu (\S+) before (\S+) .* Zsb (\S+) Zsa',t); jp=re.search(r'PAIR j projection jcl (\S+) kept \S+ of (\S+)',t)
    if not (d and pe and lg and jp): continue
    dm,E0,E1,A0,A1,aL,a1,a2,K=[float(x) for x in d.groups()]; Z,imp,gap=[float(x) for x in pe.groups()]
    Hu_e,before,Zsb=float(lg.group(2)),float(lg.group(3)),float(lg.group(4)); jcl,tot_e=float(jp.group(1)),float(jp.group(2))
    Bside=t.split('| B ')[1].split('| A ')[0]
    ax=v3(Bside,'ax'); m1=v3(Bside,'m1'); m2=v3(Bside,'m2'); e1,e2,b=frame(v3(t,'BHAT'))
    inF=lambda v: np.array([v@e1,v@e2,v@b])
    H=E0*P1+E1*P3+A0*P1@JI+A1*P3@JI+2*a1*S1[2]+2*a2*S2[2]+aL*Lv[2]; H=np.real(H+H.conj().T)/2
    w,V=np.linalg.eigh(H)
    lc=coherentL(inF(ax)); c1=coherentHalf(inF(m1*(-1)/np.linalg.norm(m1))); c2=coherentHalf(inF(m2/np.linalg.norm(m2)))
    singletPrep='PAIR prepared singlet 1' in t
    psi=np.kron(lc,singlet) if singletPrep else np.kron(lc,np.kron(c1,c2))
    o=V.T@psi; iL=1-int(round(dm)); chi=np.sum(V[4*iL:4*iL+4,:]**2,axis=0)
    B=np.abs(o)**2*chi; tot=B.sum(); Hu=np.sum(np.abs(o)**2*w)
    val.append((abs(Hu/Hu_e-1) if Hu_e!=0 else abs(Hu),abs(tot/tot_e-1)))
    B=B/tot
    # adiabatic j labels (668)
    Mk=np.rint(np.sum(V**2*Mj[:,None],axis=0)).astype(int); jl=np.zeros(12)
    for k in range(12):
        rank=sum(1 for kk in range(12) if Mk[kk]==Mk[k] and (w[kk]<w[k] or (w[kk]==w[k] and kk<k))); nh=2 if Mk[k]==0 else (1 if abs(Mk[k])==1 else 0)
        jl[k]=0.5 if rank<nh else 1.5
    RJ=B*(jl==jcl); RJ=RJ/RJ.sum() if RJ.sum()>0 else B
    Ft=before+Zsb+Z-gap
    x=(w-Ft)/u; allowed=B>1e-14
    outside=Ft<=w[allowed].min() or Ft>=w[allowed].max()
    gapout=max(w[allowed].min()-Ft,Ft-w[allowed].max(),0.0)/u
    if Ft<=w[allowed].min(): RKL=np.where(w==w[allowed].min(),1.0,0.0)
    elif Ft>=w[allowed].max(): RKL=np.where(w==w[allowed].max(),1.0,0.0)
    else:
        def mean(l):
            z=np.where(allowed,np.log(np.maximum(B,1e-300))+l*x,-np.inf); z-=z.max(); p=np.exp(z); p/=p.sum(); return p@x,p
        lo,hi=-50.0,50.0
        for _ in range(200):
            mid=0.5*(lo+hi); m,_=mean(mid)
            if m<0: lo=mid
            else: hi=mid
        _,RKL=mean(0.5*(lo+hi))
    H1=4*K*SS+2*a1*s1v[2]+2*a2*s2v[2]; w1,V1=np.linalg.eigh(np.real(H1))
    sp=V[4*iL:4*iL+4,:]                              # spin part of the m_L = dm component, unnormalized
    popk=np.array([np.abs(V1.T@(sp[:,k]/max(np.linalg.norm(sp[:,k]),1e-300)))**2 for k in range(12)])
    rho=np.zeros((4,4),complex)
    for k in range(12):
        for kk in range(12):
            nu=(w[k]-w[kk])/hP/1e6; rho+=o[k]*np.conj(o[kk])*G/(G+1j*nu)*np.outer(sp[:,k],sp[:,kk])
    rho/=np.trace(rho).real; qm=np.real(np.einsum('ai,ab,bi->i',V1,rho,V1))
    dd=int(np.argmax(np.abs(V1[3,:])**2))           # the 1S state with most |down down>
    order=np.argsort(w1)
    pk=lambda W: W@popk
    sprd=math.sqrt(max(0.0,RKL@(x*x)-(RKL@x)**2))
    res[ph].append(dict(qm=qm[order],R661=pk(B)[order],RJ=pk(RJ)[order],RKL=pk(RKL)[order],dd=int(np.where(order==dd)[0][0]),
                        kl_mean=abs(RKL@x),kl_spread=sprd,out=outside,gapout=gapout))
vmax=max(max(a,b) for a,b in val)
print(f'H669v: reconstruction of Hu and total Born weight, max relative deviation {vmax:.2e} over {len(val)} photons: {vmax<=1e-6}')
for ph,R in res.items():
    if not R: continue
    n=len(R); print(f'\n=== 1.0329 T PHEN {ph} ({"prepared singlet" if ph==1 else "product state"}): {n} photons; 1S states by energy (index of |dd>-like: {R[0]["dd"]})')
    Q=np.mean([r['qm'] for r in R],axis=0); s=np.sqrt(Q*(1-Q)/n)
    print('    QM coherent    '+' '.join(f'{x:.4f}' for x in Q)+'   2 sigma '+' '.join(f'{2*x:.4f}' for x in s))
    for rule,tag in (('R661','H669b'),('RJ','H669c'),('RKL','H669a')):
        P=np.mean([r[rule] for r in R],axis=0); ok=bool(np.all(np.abs(P-Q)<=2*s))
        print(f'    {rule:5s}          '+' '.join(f'{x:.4f}' for x in P)+f'   {tag} {ok}')
    no=[r for r in R if r['out']]; print(f'    F_t outside the spectrum of the allowed eigenstates: {len(no)}/{n}'+(f', by up to {max(r["gapout"] for r in no):.2f} ueV' if no else '')
          +f'; inside: max |<w> - F_t| {max([r["kl_mean"] for r in R if not r["out"]]+[0]):.2e} ueV')
    print(f'    RKL: max |<w> - F_t| {max(r["kl_mean"] for r in R):.2e} ueV (H669d {max(r["kl_mean"] for r in R)<=1e-3}); energy spread per photon '
          f'mean {np.mean([r["kl_spread"] for r in R]):.2f} max {max(r["kl_spread"] for r in R):.2f} ueV')
