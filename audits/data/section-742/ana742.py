#!/usr/bin/env python3
"""Audit 742 (offline, hydrogen, infinite nucleus, atomic units): engine rates OUTSIDE the EBK window (|level - n| >= 0.05) against QM.
Engine rule outside the window: s states -- counter-rotating mode at the midpoint n_m = level - s/2, L = 1, omega = s/n_m^3, weight 2 (560);
l >= 1 -- Delta l = -1 total = co-rotating modes of the initial orbit at the level (585 e_all), Delta l = +1 channels at n_m = level - s/2,
L = l + 1, weight 2(l+1)/(2l+1) (585).  Channels exist while level - s >= 2 (s states) / >= l + 2 (up).  E2 left out (<= 1e-5 of E1).
H742c: total rate / QM(n) at level n +- 0.05 (window edge; inside the window the EBK rate is QM to ~1e-4) and at n +- 0.1, 0.25 (report).
H742d: at integer levels the variant "n_m from the photon": n_m^3 = s / omega_photon, omega_photon = (1/(n-s)^2 - 1/n^2)/2 -- tau/QM against
the engine's omega_m; criterion: max |tau/QM - 1| over s, p, d, f, n = 3..8.  Run from the repository root."""
import numpy as np, os
here=os.path.dirname(os.path.abspath(__file__))
src=open(os.path.join(here,'../section-585/hyb585.py')).read().split('ns=lambda')[0]
cwd=os.getcwd(); os.chdir(os.path.join(here,'../section-585')); g={}; exec(src,g); os.chdir(cwd)
rate_qm, channels, C, AU, e_all = g['rate_qm'], g['channels'], g['C'], g['AU'], g['e_all']
f560=g['g']['fourier_components']
def cminus(nm,L,s): return f560(nm,L,s)[1]
def mid_rates(level,l,mode='m'):
    """Counter-rotating midpoint channels (s states: l = 0, L = 1, weight 2; l >= 1: up channels L = l + 1)."""
    out=[]; s=1; floor=2.0 if l==0 else l+2.0
    while level-s>=floor-1e-6:
        nm=level-0.5*s
        if mode=='g':   # n_m from the photon's omega (integer levels only)
            om=0.5*(1.0/(level-s)**2-1.0/level**2); nm=(s/om)**(1.0/3.0)
        L=1.0 if l==0 else l+1.0; w=2.0 if l==0 else 2.0*(l+1)/(2*l+1)
        om=s/nm**3
        out.append(w*(2/3)*om**3*cminus(nm,L,s)**2/C**3); s+=1
    return out
def engine_total(level,l,mode='m'):
    t=sum(mid_rates(level,l,mode))
    if l>=1: t+=e_all(level,l,N=1<<16)
    return t
def qm_total(n,l): return sum(rate_qm(n,l,n2,l2) for n2,l2 in channels(n,l))
names="spdfgh"
print("H742c: total engine rate outside the window / QM(n)   (0 = no channel open at that level)")
print("state   tau_QM[ns]  " + "  ".join(f"{d:+6.2f}" for d in (-0.25,-0.1,-0.05,0.05,0.1,0.25)) + "   rule@n (window off)")
for n in range(2,7):
    for l in range(0,n):
        if l==0 and n<3: continue
        q=qm_total(n,l)
        if q<=0: continue
        r=[engine_total(n+d,l)/q for d in (-0.25,-0.1,-0.05,0.05,0.1,0.25)]
        print(f"{n}{names[l]}    {1e9/(q*AU):9.3f}   " + "  ".join(f"{x:6.3f}" for x in r) + f"   {engine_total(float(n),l)/q:6.3f}")
print("\nH742d: tau/QM at integer levels, omega_m (engine outside the window) vs n_m from the photon's omega")
worst={'m':0,'g':0}
print("state   omega_m   photon-n_m")
for n in range(3,9):
    for l in range(0,4):
        if l>=n: continue
        q=qm_total(n,l)
        if l>=1 and n-1<l+2 and l>=1:   # no up channel: variants identical
            pass
        t={k:q/engine_total(float(n),l,k) for k in ('m','g')}
        for k in t: worst[k]=max(worst[k],abs(t[k]-1))
        print(f"{n}{names[l]}     {t['m']:6.3f}    {t['g']:6.3f}")
print(f"max |tau/QM - 1|: omega_m {worst['m']:.3f}, photon-n_m {worst['g']:.3f}")
