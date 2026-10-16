# Audit 485: a carrier for the fast cone phase in 1e13-1e16 rad/s, or ZPF sustaining the transverse amplitude without net radiation.
import math
hbar=1.054571817e-34; c=299792458.0; me=9.1093837015e-31; e=1.602176634e-19; mu0=1.25663706212e-6; muB=9.2740100783e-24; eps0=8.8541878128e-12
g=2.00231930436; gam=g*e/(2*me); re=e*e/(4*math.pi*eps0*me*c*c); a0=5.29177210903e-11
mperp=math.sqrt(2)*muB; t=hbar/math.sqrt(2)
life={'p-Ps 1s':125.14e-12,'Ps 2p':3.19e-9,'o-Ps 1s':142.04e-9}
print('(1) radiation and (2) drive, per omega_f:')
print('  omega_f [rad/s]  P [W]        tau_rad [s]   tau_rad/tau(o-Ps)  B_f = omega/gamma [T]')
for w in (1e13,1e14,1e15,1e16):
    P=mu0*w**4*mperp**2/(6*math.pi*c**3); E=hbar*w*(t/hbar); tau=E/P
    print(f'  {w:.0e}            {P:.3e}   {tau:.3e}     {tau/life["o-Ps 1s"]:.3e}          {w/gam:.3e}')
fields={'contact Ps 1s (3.47e11 rad/s, 479/470)':3.4667e11,'spin-orbit Ps 1s classical (1.654e12, 477)':1.654e12,'spin-orbit Ps 2p (2.297e10, 467)':2.297e10,
        'spin-orbit H 1s classical (4.41e12, 469)':4.41e12}
print('   fields available in the atom (precession rate / gamma):')
for k,v in fields.items(): print(f'     {k:46s} {v/gam:8.2f} T  ({v:.2e} rad/s)')
print('\n(3) reservoir (independent of omega): hidden angular momentum (r_e/2r) hbar vs sqrt2 hbar -- Ps 1s 9.4e-6, Ps 2p 2.4e-6, H 1s 1.9e-5 (481)')
print('\n(4) ZPF (SED) equilibrium of a precession mode: exponential cap action, <A> = hbar/2, single-sided; approach time ~ 1/Gamma_rad:')
for name,w in (('contact Ps 1s',3.4667e11),('spin-orbit H 1s',4.41e12),('omega_f 1e13',1e13),('omega_f 1e16',1e16)):
    P=mu0*w**4*mperp**2/(6*math.pi*c**3); G=P/(hbar*w*(t/hbar))
    print(f'   {name:16s} omega {w:.2e}: Gamma_rad {G:.3e} 1/s, approach time {1/G:.3e} s = {1/G/life["o-Ps 1s"]:.2e} x tau(o-Ps)')
print('   distribution: P(A) = (2/hbar) exp(-2A/hbar) -- with |S| = hbar (Langer) S_z = hbar - A: one broad line on the field side, mean hbar/2;')
print('   anti-aligned orientations relax radiatively -- no second line; fixed delta is not a stationary state, only its mean is.')
print('\n(5) orbital rotation as the carrier:')
for name,n,r in (('Ps 1s',1,2*a0),('Ps 2p',2,8*a0),('H 1s',1,a0)):
    mu=me/2 if 'Ps' in name else me*1836.15267/(1836.15267+1); k=e*e/(4*math.pi*eps0)
    worb=math.sqrt(k/(mu*r**3)); ratio=0.75*re/r
    print(f'   {name}: omega_orb {worb:.3e} rad/s (in the window); Omega_SO/omega_orb = (3/4) r_e/r = {ratio:.2e}; torque needed to co-rotate the spin '
          f'{worb*t:.2e} J vs available {ratio*worb*t:.2e} J (factor {1/ratio:.1e})')
