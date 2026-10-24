# Audit 619 step 1: Born probabilities of coherent states on coupled states vs rule-Q basins; isotropic averages.
import numpy as np, math
from math import factorial as fa
def cg(j1,m1,j2,m2,J,M):
    if abs(m1+m2-M)>1e-9 or J<abs(j1-j2)-1e-9 or J>j1+j2+1e-9 or abs(M)>J+1e-9: return 0.0
    f=lambda x: fa(int(round(x)))
    pre=math.sqrt((2*J+1)*f(J+j1-j2)*f(J-j1+j2)*f(j1+j2-J)/f(j1+j2+J+1))*math.sqrt(f(J+M)*f(J-M)*f(j1-m1)*f(j1+m1)*f(j2-m2)*f(j2+m2))
    s=0.0
    for k in range(0,100):
        a=[k, j1+j2-J-k, j1-m1-k, j2+m2-k, J-j2+m1+k, J-j1-m2+k]
        if any(x< -1e-9 for x in a): continue
        if a[1]<-1e-9 or a[2]<-1e-9 or a[3]<-1e-9: break
        s+=(-1)**k/np.prod([f(x) for x in a])
    return pre*s
def dprob(j,m,theta):   # |d^j_{m,j}(theta)|^2 : coherent state along theta (max projection j) on |j m> about the reference axis
    c=math.cos(theta/2); s=math.sin(theta/2)
    return math.comb(int(2*j),int(round(j-m)))*c**(2*(j+m))*s**(2*(j-m))
def pborn(j1,j2,theta):
    # first moment coherent along the axis (m1 = j1), second along angle theta to it
    out={}
    for J in np.arange(abs(j1-j2),j1+j2+0.1,1.0):
        p=0.0
        for m2 in np.arange(-j2,j2+0.1,1.0): p+=dprob(j2,m2,theta)*cg(j1,j1,j2,m2,J,j1+m2)**2
        out[round(J,1)]=p
    return out
def pbasin(j1,j2,theta,len1,len2):
    tot=math.sqrt(len1**2+len2**2+2*len1*len2*math.cos(theta))
    J=round(tot-0.5) if (j1+j2)%1==0 else round(tot-0.5-0.5)+0.5
    J=min(max(J,abs(j1-j2)),j1+j2); return {round(J,1):1.0}
th=np.linspace(0,math.pi,20001); wt=np.sin(th); wt/=np.trapezoid(wt,th)
worst=0
for tag,(j2,len2) in (("j (H, spin 1/2)",(0.5,1.0)),("J (Ps triplet, S = 1)",(1.0,1.5))):
    for l in range(1,5):
        Ps=[pborn(l,j2,t) for t in th]; Js=sorted(Ps[0])
        iso={J:np.trapezoid([p[J] for p in Ps]*wt,th) for J in Js}
        ref={J:(2*J+1)/((2*j2+1)*(2*l+1)) for J in Js}
        d=max(abs(iso[J]-ref[J]) for J in Js); worst=max(worst,d)
        Pb=[pbasin(l,j2,t,l+0.5,len2) for t in th]
        diff=np.trapezoid([sum(abs(p.get(J,0)-q.get(J,0)) for J in Js)/2 for p,q in zip(Ps,Pb)]*wt,th)
        isob={J:np.trapezoid([q.get(J,0) for q in Pb]*wt,th) for J in Js}
        print(f"{tag} l={l}: Born isotropic "+" ".join(f"{J}:{iso[J]:.6f}" for J in Js)+f"  exact (2J+1)/... max dev {d:.1e}; basin isotropic "+" ".join(f"{J}:{isob[J]:.4f}" for J in Js)+f"; per-configuration mean TV(Born, basin) {diff:.3f}")
# F: j coherent (max projection) x nuclear spin 1/2 -- same as spin-1/2 coupling with j1 = j
for j in (0.5,1.5,2.5):
    Ps=[pborn(j,0.5,t) for t in th]; Fs=sorted(Ps[0]); iso={F:np.trapezoid([p[F] for p in Ps]*wt,th) for F in Fs}
    d=max(abs(iso[F]-(2*F+1)/(2*(2*j+1))) for F in Fs); worst=max(worst,d)
    print(f"F (j={j}, I=1/2): Born isotropic "+" ".join(f"{F}:{iso[F]:.6f}" for F in Fs)+f"  max dev {d:.1e}")
# Ps 1S in a field: isotropic triplet rho = P_T/3 projected on eigenstates -> T+1, T-1 one third each, block M = 0 one third (eigenvectors within the block)
print(f"1S triplet isotropic: P(M=+1) = P(M=-1) = P(M=0 block) = 1/3 exactly (P_T/3 is diagonal in M with 1/3 per triplet state)")
print(f"H619a max deviation {worst:.1e}")
