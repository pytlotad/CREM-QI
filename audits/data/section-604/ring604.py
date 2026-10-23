# Audit 604: internal frequency shift of a classical charged ring (electron as a loop) moving on a Kepler orbit in the proton's Coulomb field.
# Vectorized over orientations (6), internal speeds (3) and eccentricities (2); hbar = m = c = 1, k = alpha, R = 1.
import numpy as np, sys
alpha=7.2973525693e-3
orients=np.array([(0,0,1),(1,0,0),(0,1,0),(1,1,1),(1,-1,0.3),(0.2,0.7,-1)],float)
vs=np.array([0.3,0.6,0.9]); es=np.array([0.0,0.5])
def run(n):
    R=1.0; a=n*n/alpha; Om=np.sqrt(alpha/a**3); T=2*np.pi/Om
    combos=[(ie,iv,io) for ie in range(2) for iv in range(3) for io in range(6)]
    E=np.array([es[c[0]] for c in combos]); V=np.array([vs[c[1]] for c in combos]); NH=orients[[c[2] for c in combos]]
    NH/=np.linalg.norm(NH,axis=1,keepdims=True)
    T1=np.where(np.abs(NH[:,:1])<0.9,np.array([[1.,0,0]]),np.array([[0,1.,0]]))
    E1=np.cross(NH,T1); E1/=np.linalg.norm(E1,axis=1,keepdims=True); E2=np.cross(NH,E1)
    P0=R*V/np.sqrt(1-V*V); W0=V/R
    dt=2*np.pi/vs.max()/20; N=int(np.ceil(T/dt)); dt=T/N
    ustate={}
    def rc(t,key):
        M=Om*t; u=ustate.get(key,np.full(2,M))
        for _ in range(3): u=u-(u-es*np.sin(u)-M)/(1-es*np.cos(u))
        ustate[key]=u
        pos=np.stack([a*(np.cos(u)-es), a*np.sqrt(1-es*es)*np.sin(u), np.zeros(2)],1)   # (2,3) per eccentricity
        return pos[[c[0] for c in combos]]
    def f(th,p,t,key):
        c=rc(t,key); cth=np.cos(th)[:,None]; sth=np.sin(th)[:,None]
        rho=R*(cth*E1+sth*E2); drho=R*(-sth*E1+cth*E2); x=c+rho
        Ex=alpha*x/np.linalg.norm(x,axis=1,keepdims=True)**3; Ec=alpha*c/np.linalg.norm(c,axis=1,keepdims=True)**3
        dU=np.sum((Ex-Ec)*drho,1)
        return (p/R**2)/np.sqrt(1+p*p/R**2), -dU
    th=np.zeros(len(combos)); p=P0.copy(); t=0.0
    for i in range(N):
        a1,b1=f(th,p,t,'0'); a2,b2=f(th+0.5*dt*a1,p+0.5*dt*b1,t+0.5*dt,'h'); a3,b3=f(th+0.5*dt*a2,p+0.5*dt*b2,t+0.5*dt,'h'); a4,b4=f(th+dt*a3,p+dt*b3,t+dt,'1')
        th=th+dt/6*(a1+2*a2+2*a3+a4); p=p+dt/6*(b1+2*b2+2*b3+b4); t+=dt
        ustate['0']=ustate['1']
    sh=th/(W0*T)-1.0; req=alpha**2/n**2
    for ie in range(2):
        for iv in range(3):
            idx=[k for k,c in enumerate(combos) if c[0]==ie and c[1]==iv]
            m=np.max(np.abs(sh[idx]))
            print(f"{n}  {es[ie]:.1f}  {vs[iv]:.1f}   max|<dw/w>| {m:.3e} (mean {np.mean(sh[idx]):+.3e})   required {-req:.3e}   ratio {m/req:.1e}",flush=True)
print("n  e    v     classical ring shift over one orbit vs Bohr requirement <V>/mc^2")
for n in (1,2): run(n)
print(f"scale alpha (R/a)^3 at n = 1: {alpha**4:.2e};  squared: {alpha**8:.2e};  required alpha^2: {alpha**2:.2e}")
