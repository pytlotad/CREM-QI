# Audit 603: the factor 4 of the internal clock (audit 380) -- proper-time slip vs Lagrangian phase, per particle and summed.
import numpy as np
alpha=7.2973525693e-3
def orbit(e,n,mu,k,N=1<<16):
    # Kepler relative orbit with hbar = 1, a = n^2 hbar^2/(mu k), J = n (principal action)
    a=n*n/(mu*k); Om=np.sqrt(k/(mu*a**3)); M=np.arange(N)*2*np.pi/N; u=np.full_like(M,np.pi)
    for _ in range(200): u-=(u-e*np.sin(u)-M)/(1-e*np.cos(u))
    r=a*(1-e*np.cos(u)); udot=Om/(1-e*np.cos(u))
    v2=(a*udot)**2*(np.sin(u)**2+(1-e*e)*np.cos(u)**2)
    dt=(2*np.pi/Om)/N
    return r,v2,dt,a,Om
for name,m1,m2 in (("Ps",1.0,1.0),("H",1.0,1836.15267343)):
    mu=m1*m2/(m1+m2); c=1/alpha; k=1.0              # atomic-like units: hbar = 1, e^2 = k = 1, c = 1/alpha
    for n in (1,2,3):
        for e in (0.0,0.5,0.9):
            r,v2,dt,a,Om=orbit(e,n,mu,k); J=n
            vi2=[v2*(m2/(m1+m2))**2, v2*(m1/(m1+m2))**2]   # particle speeds in the CM frame
            slips=[m*c*c*np.sum(1-np.sqrt(1-w/c**2))*dt/(2*np.pi) for m,w in zip((m1,m2),vi2)]   # cycles per orbit
            # (B) nonrelativistic Lagrangian phases summed vs reference E: sum_i oint (T_i) dt - oint V dt + E T = oint p.dq
            T=0.5*mu*v2; V=-k/r; E=-k/(2*a); P=2*np.pi/Om
            B_=(np.sum(T-V)*dt+E*P)/(2*np.pi*J)
            print(f"{name} n={n} e={e}: proper-time slip / J: clock 1 {slips[0]/J:.6f}, clock 2 {slips[1]/J:.6f};  (B) Lagrangian sum vs E: {B_:.12f}")
