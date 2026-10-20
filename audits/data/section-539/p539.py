# Audit 539: Ps 1S m = 0 in a field, two-spin model.  Classical Zeeman U_Z = Z(n1z - n2z) vs rule U_Z = Z(n1z - n2z)|n1 + n2|/2; contact K n1.n2, K = dE/4.
import numpy as np
h=6.62607015e-34; hbar=h/(2*np.pi); muB=9.2740100783e-24; g=2.00231930436092; dE=203.38910e9*h; K=dE/4
def rates(n1,n2,Z,rule):
    if not rule:
        return (2/hbar)*(K*n2+np.array([0,0,Z])), (2/hbar)*(K*n1-np.array([0,0,Z]))
    S=n1+n2; s=np.linalg.norm(S)/2; Dz=n1[2]-n2[2]
    common=Z*Dz*S/(4*max(s,1e-300))
    return (2/hbar)*(K*n2+np.array([0,0,Z*s])+common), (2/hbar)*(K*n1-np.array([0,0,Z*s])+common)
def run(B,start,rule,T=1e-9,dt=2.5e-14):
    Z=(g/2)*muB*B; x=4*Z/dE
    if start=='parallel': n1=np.array([1.0,0,0]); n2=n1.copy()
    else:
        phi=0.5*np.arctan(x) if rule else np.arcsin(min(1.0,Z/(2*K)))
        n1=np.array([np.cos(phi),0,np.sin(phi)]); n2=np.array([np.cos(phi),0,-np.sin(phi)])
    def f(a,b):
        o1,o2=rates(a,b,Z,rule); return np.cross(o1,a),np.cross(o2,b)
    wsum=0.0; N=int(T/dt)
    for i in range(N):
        a1,a2=f(n1,n2); b1,b2=f(n1+0.5*dt*a1,n2+0.5*dt*a2); c1,c2=f(n1+0.5*dt*b1,n2+0.5*dt*b2); d1,d2=f(n1+dt*c1,n2+dt*c2)
        n1=n1+dt/6*(a1+2*b1+2*c1+d1); n2=n2+dt/6*(a2+2*b2+2*c2+d2); n1/=np.linalg.norm(n1); n2/=np.linalg.norm(n2)
        wsum+=(1-n1@n2)/2
    return wsum/N
for B in (0.1,0.3,1.0,3.0,10.0):
    x=2*g*muB*B/dE; PS=(1-1/np.sqrt(1+x*x))/2; T0avg=x*x/(2*(1+x*x))
    cs=run(B,'stationary',False); rs=run(B,'stationary',True); cp=run(B,'parallel',False); rp=run(B,'parallel',True)
    print(f"B {B:4.1f} T x {x:.4f}: P_S {PS:.5e} | (A) stationary: classical {cs/PS:.4f}, rule {rs/PS:.4f} x P_S | (B) parallel start: QM T0 average {T0avg:.5e};"
          f" classical {cp:.5e} ({cp/T0avg:.4f}), rule {rp:.5e} ({rp/T0avg:.4f})",flush=True)
