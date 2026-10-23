# Audit 583: s-state E1 (560) and E2 (583) mode-rule rates of hydrogen in atomic units (infinite-mass approximation as in the 560 script).
import numpy as np, math
C=137.035999; AU=4.134137333e16
def comp(n_orb,L,s,quad,N=8192):
    a=n_orb**2; e=math.sqrt(max(0,1-(L/n_orb)**2)); b=a*math.sqrt(1-e*e)
    u=(np.arange(N)+0.5)*2*np.pi/N; M=u-e*np.sin(u); z=a*(np.cos(u)-e)+1j*b*np.sin(u); w=(1-e*np.cos(u))/N
    f=z*z if quad else z
    return abs(np.sum(f*np.exp(1j*s*M)*w))
for n in range(3,7):
    e1=sum(2*(2/3)*(s/((n-0.5*s)**3))**3/C**3*comp(n-0.5*s,1.0,s,False)**2 for s in range(1,n-1))
    e2=0.0
    for s in range(1,n):
        if n-s<3: continue
        nm=n-0.5*s; om=s/nm**3
        e2+=4*om**6*comp(nm,1.5,s,True)**2/(40*C**5)/om
    print(f"{n}s: E1 rate {e1*AU:.4e} /s (tau {1e9/(e1*AU):.2f} ns), E2 rule {e2*AU:.4e} /s, E2/E1 {e2/e1:.3e}")
