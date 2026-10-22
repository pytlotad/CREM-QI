# Audit 556 H3/H4: (H3) the classical Fourier spectrum of the hyperbolic passage per L, photons dE/domega/omega assigned to n; (H4) the stationary-phase rule with phi' = L/r^2.
import math, numpy as np
exec(open('../section-555/pw555.py').read().split("for KoverB,label in")[0])
def fourier_P(L,K,nmax):
    ah=1/(2*K); e=math.sqrt(1+2*K*L*L); T=math.sqrt(ah**3)
    F=np.linspace(-25,25,400001); t=T*(e*np.sinh(F)-F); x=ah*(e-np.cosh(F)); y=ah*math.sqrt(e*e-1)*np.sinh(F)
    r=np.sqrt(x*x+y*y); ax=-x/r**3; ay=-y/r**3; dtdF=T*(e*np.cosh(F)-1)
    out=np.zeros(nmax+1)
    for n in range(1,nmax+1):
        lo=K+0.5/((n+0.5)**2); hi=K+0.5/((n-0.5)**2) if n>1 else K+2.0   # band of n: n_cl in [n-1/2, n+1/2)
        ws=np.linspace(lo,hi,60); acc=0.0
        for i,w in enumerate(ws):
            ph=np.exp(1j*w*t)*dtdF; Ax=np.trapezoid(ax*ph,F); Ay=np.trapezoid(ay*ph,F)
            dEdw=(2*alpha**3/(3*math.pi))*(abs(Ax)**2+abs(Ay)**2)
            acc+=dEdw/w
        out[n]=acc*(ws[1]-ws[0])
    return out
def phidot_P(L,K,nmax):
    ah=1/(2*K); e=math.sqrt(1+2*K*L*L); F=np.linspace(-12,12,60001); u=e*np.cosh(F)-1; r=ah*u
    om=L/r**2; P=(2/3)*alpha**3*r**-4; dN=P/om*np.sqrt(ah**3)*u; dF=F[1]-F[0]
    out=np.zeros(nmax+1)
    for i in range(len(F)):
        if not om[i]>K: continue
        ncl=math.sqrt(0.5/(om[i]-K))
        if ncl<0.5: continue
        n=max(1,int(round(ncl)))
        if n<=nmax: out[n]+=dN[i]*dF
    return out
KoverB=0.18; K=KoverB/2; kw=math.sqrt(2*K); nmax=4
print(f"Ps K/B = {KoverB}: model/QM per partial wave, n = 1..4 (band average over L, weight L dL, 9 points)")
tot={'QM':0,'H3':0,'H4':0,'549':0}
for lin in range(0,3):
    q=np.array([qm_P(lin,n,kw) for n in range(1,nmax+1)])
    Ls=np.linspace(lin+0.06,lin+0.94,9)
    h3=sum(fourier_P(L,K,nmax)*L for L in Ls)[1:]/Ls.sum()
    h4=sum(phidot_P(L,K,nmax)*L for L in Ls)[1:]/Ls.sum()
    r549=sum(cl_P(L,K,nmax,False)*L for L in Ls)[1:]/Ls.sum()
    for name,v in (('549 (circular omega)',r549),('H3 Fourier',h3),('H4 phi-dot',h4)):
        print(f"  l_in={lin} {name:22s}: "+' '.join((f"{a/b:.2f}" if b>0 else (f"[{a:.1e}]" if a>0 else '-')) for a,b in zip(v,q)))
    w=2*lin+1; tot['QM']+=w*q.sum(); tot['H3']+=w*h3.sum(); tot['H4']+=w*h4.sum(); tot['549']+=w*r549.sum()
print("  sum over l_in = 0-2, n = 1-4 (weights 2l+1): model/QM -- 549 %.3f, H3 %.3f, H4 %.3f"%(tot['549']/tot['QM'],tot['H3']/tot['QM'],tot['H4']/tot['QM']))
