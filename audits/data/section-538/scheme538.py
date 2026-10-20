# Audit 538: the engine's s-band spin scheme in field on the two-spin model (Ps 1S, magnetic contact K_c, B = 1 T): explicit midpoint for the internal
# (contact) rates with exact rotations, then an exact external rotation (503 splitting); step bounded by an angle per substep.  Variants to find the error.
import numpy as np, sys
h=6.62607015e-34; hbar=h/(2*np.pi); muB=9.2740100783e-24; g=2.00231930436092
Kc=204.6577e9*h/1.748264/4
def rot(v,w,dt):
    th=np.linalg.norm(w)*dt
    if th==0: return v
    k=w/np.linalg.norm(w); return v*np.cos(th)+np.cross(k,v)*np.sin(th)+k*(k@v)*(1-np.cos(th))
def rates(n1,n2,K,Z):   # rate vectors (dn/dt = Omega x n), internal and external parts
    i1=(2/hbar)*K*n2; i2=(2/hbar)*K*n1; e1=(2/hbar)*np.array([0,0,Z]); e2=-e1
    return i1,i2,e1,e2
def run(scheme,angle,B=1.0,T=1e-9,K=Kc):
    Z=(g/2)*muB*B; n1=np.array([1.0,0,0]); n2=n1.copy(); t=0.0; wsum=0.0; tw=0.0
    while t<T:
        i1,i2,e1,e2=rates(n1,n2,K,Z)
        sp=max(np.linalg.norm(i1+e1),np.linalg.norm(i2+e2)); dt=min(angle/sp,T-t)
        if scheme=='engine':        # midpoint internal with exact rotations, then exact external by the midpoint external rate
            m1=rot(n1,i1,dt/2); m2=rot(n2,i2,dt/2); j1,j2,f1,f2=rates(m1,m2,K,Z)
            n1=rot(n1,j1,dt); n2=rot(n2,j2,dt); n1=rot(n1,f1,dt); n2=rot(n2,f2,dt)
        elif scheme=='strang':      # external half, internal midpoint, external half
            n1=rot(n1,e1,dt/2); n2=rot(n2,e2,dt/2)
            i1,i2,_,_=rates(n1,n2,K,Z); m1=rot(n1,i1,dt/2); m2=rot(n2,i2,dt/2); j1,j2,_,_=rates(m1,m2,K,Z)
            n1=rot(n1,j1,dt); n2=rot(n2,j2,dt); n1=rot(n1,e1,dt/2); n2=rot(n2,e2,dt/2)
        elif scheme=='unsplit':     # midpoint with the total rate
            m1=rot(n1,i1+e1,dt/2); m2=rot(n2,i2+e2,dt/2); j1,j2,f1,f2=rates(m1,m2,K,Z)
            n1=rot(n1,j1+f1,dt); n2=rot(n2,j2+f2,dt)
        n1/=np.linalg.norm(n1); n2/=np.linalg.norm(n2)
        wsum+=(1-n1@n2)/2*dt; tw+=dt; t+=dt
    return wsum/tw
for scheme in ('engine','strang','unsplit'):
    print(scheme, ' '.join(f"angle {a}: <w> {run(scheme,a):.4e}" for a in (0.05,0.01,0.002)), flush=True)
print("reference (RK4, 536/537): 0.1185")
