# Audit 536: classical two-spin model of Ps 1S in a field vs the engine (535) and QM mixing; variants of the coupling and of the starting configuration.
import numpy as np
h=6.62607015e-34; hbar=h/(2*np.pi); muB=9.2740100783e-24; g=2.00231930436092
dE=203.38910e9*h                      # measured (Ritter et al. 1984)
dE_model=204.6577e9*h                 # engine's 1S HFS (p535)
Kc=dE_model*(4/7)/4                   # engine torque: magnetic contact only (annihilation 3/7 has no torque)
Kr=dE/4                               # rule B: slope of the level energy at the midpoint
def run(K,B,start,T=1e-9,dt=2.5e-14):
    Z=(g/2)*muB*B
    if start=='parallel': n1=np.array([1.0,0,0]); n2=n1.copy()
    else:
        s=min(1.0,Z/(2*K)); c=np.sqrt(1-s*s); n1=np.array([c,0,s]); n2=np.array([c,0,-s])
    def f(n1,n2):
        g1=(2/hbar)*(K*n2+np.array([0,0,Z])); g2=(2/hbar)*(K*n1-np.array([0,0,Z]))
        return np.cross(g1,n1),np.cross(g2,n2)
    N=int(T/dt); wsum=0.0
    for i in range(N):
        a1,a2=f(n1,n2); b1,b2=f(n1+0.5*dt*a1,n2+0.5*dt*a2); c1,c2=f(n1+0.5*dt*b1,n2+0.5*dt*b2); d1,d2=f(n1+dt*c1,n2+dt*c2)
        n1=n1+dt/6*(a1+2*b1+2*c1+d1); n2=n2+dt/6*(a2+2*b2+2*c2+d2); n1/=np.linalg.norm(n1); n2/=np.linalg.norm(n2)
        wsum+=(1-n1@n2)/2
    return wsum/N
print("check (i): relative mode at B = 0 precesses at 4K/hbar -> K_c gives %.3f GHz, K = dE/4 gives %.3f GHz"%(4*Kc/h*1e-9,4*Kr/h*1e-9))
for B in (0.1,0.3,1.0,3.0):
    x=2*g*muB*B/dE; y=x/(1+np.sqrt(1+x*x)); PS=y*y/(1+y*y)
    A=run(Kc,B,'parallel'); Bv=run(Kr,B,'parallel'); C=run(Kr,B,'stationary'); D=run(Kc,B,'stationary')
    print(f"B {B:.1f} T: QM P_S {PS:.4e} | (A) K_c, parallel {A:.4e} ({A/PS:.3f}) | (B) dE/4, parallel {Bv:.4e} ({Bv/PS:.3f}) | (C) dE/4, stationary {C:.4e} ({C/PS:.3f})"
          f" | (D) K_c, stationary {D:.4e} ({D/PS:.3f}) | classical stationary w = (Z/2K)^2 with dE/4: {((g/2)*muB*B/(2*Kr))**2:.4e}")
print("engine (535): <w>/P_S = 6.06 / 6.01 / 16.3 / 5.70")
