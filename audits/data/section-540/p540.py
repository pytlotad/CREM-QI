# Audit 540: two ways to remove the factor 2 (nutation) for o-Ps in a field -- (1) start rule (rule Q basin + stationary state), (2) damping from decay -- vs exact QM {S, T0} with decay.
import numpy as np
h=6.62607015e-34; hbar=h/(2*np.pi); muB=9.2740100783e-24; g=2.00231930436092; dE=203.38910e9*h; K=dE/4; J=dE/4
Gp=1/125.2e-12; Go=Gp*4*(np.pi**2-9)*7.2973525693e-3/(9*np.pi); Gd=Gp-Go
def qm_survival(B,times,psi0):
    Z=(g/2)*muB*B; Hm=np.array([[-3*J-0.5j*hbar*Gp,2*Z],[2*Z,J-0.5j*hbar*Go]])/hbar
    w,V=np.linalg.eig(Hm); Vi=np.linalg.inv(V); c=Vi@psi0
    out=[]
    for t in times:
        psi=V@(c*np.exp(-1j*w*t)); out.append(np.vdot(psi,psi).real)
    late=min(-2*np.imag(w))          # slowest decay rate of the mixed pair
    return np.array(out),late
def rates(n1,n2,Z):
    S=n1+n2; s=np.linalg.norm(S)/2; Dz=n1[2]-n2[2]; common=Z*Dz*S/(4*max(s,1e-300))
    return (2/hbar)*(K*n2+np.array([0,0,Z*s])+common),(2/hbar)*(K*n1-np.array([0,0,Z*s])+common)
def classical(B,start,theta,times,damp,dt=2.5e-14):
    Z=(g/2)*muB*B; x=4*Z/dE
    if start=='stationary': phi=0.5*np.arctan(x); n1=np.array([np.cos(phi),0,np.sin(phi)]); n2=np.array([np.cos(phi),0,-np.sin(phi)])
    elif start=='alongB': n1=np.array([0,0,1.0]); n2=n1.copy()
    else: n1=np.array([np.sin(theta),0,np.cos(theta)]); n2=n1.copy()
    kap=Gd/4 if damp else 0.0
    def f(a,b):
        o1,o2=rates(a,b,Z); d1=np.cross(o1,a)+kap*(b-(a@b)*a); d2=np.cross(o2,b)+kap*(a-(a@b)*b); return d1,d2
    t=0.0; lnS=0.0; out=[]; k=0; tw=[]
    T=times[-1]
    while k<len(times):
        a1,a2=f(n1,n2); b1,b2=f(n1+0.5*dt*a1,n2+0.5*dt*a2); c1,c2=f(n1+0.5*dt*b1,n2+0.5*dt*b2); d1,d2=f(n1+dt*c1,n2+dt*c2)
        n1=n1+dt/6*(a1+2*b1+2*c1+d1); n2=n2+dt/6*(a2+2*b2+2*c2+d2); n1/=np.linalg.norm(n1); n2/=np.linalg.norm(n2)
        w=(1-n1@n2)/2; lnS-=(Go+w*Gd)*dt; t+=dt; tw.append(w)
        if t>=times[k]: out.append(np.exp(lnS)); k+=1
    wl=np.mean(tw[-int(0.2e-9/dt):])
    return np.array(out),Go+wl*Gd
times=[0.5e-9,2e-9]
for B in (0.1,1.0,3.0):
    x=4*(g/2)*muB*B/dE; PS=(1-1/np.sqrt(1+x*x))/2
    Sq,late=qm_survival(B,times,np.array([0,1.0+0j]))
    S1,l1=classical(B,'stationary',None,times,False); S2,l2=classical(B,'perp',np.pi/2,times,True); S0,l0=classical(B,'perp',np.pi/2,times,False)
    print(f"B {B} T: QM T0 start S(0.5, 2 ns) {Sq[0]:.6f} {Sq[1]:.6f}, late rate {late:.5e}/s (eigen: Go + P_S Gd = {Go+PS*Gd:.5e})")
    print(f"   (1) start rule (stationary): S {S1[0]:.6f} {S1[1]:.6f}, late {l1:.5e} ({l1/late:.4f}) | (2) damping, start T0-like: S {S2[0]:.6f} {S2[1]:.6f}, late {l2:.5e} ({l2/late:.4f})"
          f" | no fix: S {S0[0]:.6f} {S0[1]:.6f}, late {l0:.5e} ({l0/late:.4f})",flush=True)
# isotropic ensemble at 1 T: QM late quench (excess over Go) = (1/3) P_S Gd ; (1) start rule by construction 1/3 x stationary; (2) damping vs angle
B=1.0; x=4*(g/2)*muB*B/dE; PS=(1-1/np.sqrt(1+x*x))/2
print(f"ensemble 1 T: QM excess late rate (1/3) P_S Gd = {PS*Gd/3:.5e}/s; (1) start rule: {PS*Gd/3:.5e} (m = 0 basin 1/3, stationary w = P_S)")
ex=[]
for th in (15,30,45,60,75,90):
    _,l=classical(B,'perp',np.radians(th),[1.0e-9],True); ex.append((th,l-Go)); print(f"   (2) damping, parallel spins at {th} deg to B: late excess {l-Go:.5e} (/(P_S Gd) = {(l-Go)/(PS*Gd):.4f}; sin^2 = {np.sin(np.radians(th))**2:.4f})",flush=True)
