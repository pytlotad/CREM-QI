# Audit 481: phase coupling of the transverse cone components of two spins without import 8 and rule Q.
import math, random
hbar=1.054571817e-34; c=299792458.0; me=9.1093837015e-31; e=1.602176634e-19; eps0=8.8541878128e-12; muB=9.2740100783e-24
a0=5.29177210903e-11; re=e*e/(4*math.pi*eps0*me*c*c); kC=e*e/(4*math.pi*eps0)
t=1/math.sqrt(2)   # transverse spin component in hbar
print('(1) reservoir: hidden angular momentum available (r_e/(2r)) hbar vs required swing sqrt(2) hbar per fast period')
for name,r in (('Ps 2p (4 a_Ps)',8*a0),('Ps 1s (a_Ps)',2*a0),('H 1s (a0)',a0),('r* = 193 fm',193.3e-15),('r_e',re),('0.7 r_e',0.7*re)):
    avail=muB*kC/(e*r*c*c*hbar)*e   # |m||E| r / c^2 = mu_B (kC e/(e r^2)) r / c^2 ... in hbar
    avail=muB*(kC/(e*r*r))*r/(c*c)/hbar
    print(f'   {name:16s} r = {r:.3e} m: available {avail:.3e} hbar = {avail/(2*t):.3e} of the required {2*t:.3f} hbar')
# (2) conservative coupling: H = J [S1z S2z + t1 t2 cos(psi1 - psi2)], t_i = sqrt(3/4 - S_iz^2); (psi_i, S_iz) canonical pairs
print('\n(2) conservative (Heisenberg-like) coupling in canonical variables (psi_i, S_iz):')
def Hrel(x,p):   # relative phase x = psi1 - psi2, with S1z = -S2z = p (opposite projections, total S_z = 0 conserved)
    t1=math.sqrt(max(0,0.75-p*p)); return -p*p+t1*t1*math.cos(x)
# Jacobian of Hamilton's equations at fixed points (x, p) = (0, 0) and (pi, 0): dx/dt = dH/dp, dp/dt = -dH/dx (with factor 2 from two spins absorbed)
for x0 in (0.0,math.pi):
    hpp=(Hrel(x0,1e-4)-2*Hrel(x0,0)+Hrel(x0,-1e-4))/1e-8; hxx=(Hrel(x0+1e-4,0)-2*Hrel(x0,0)+Hrel(x0-1e-4,0))/1e-8
    kind='centre (oscillation, no attractor)' if hpp*hxx>0 else 'saddle (unstable)'
    print(f'   dpsi = {x0:.3f}: H_pp {hpp:+.3f}, H_xx {hxx:+.3f} -> {kind}; trace of the linearized flow = 0 (Hamiltonian, area preserving)')
# (3) coherent magnetic radiation reaction in the near zone on the rotating transverse moments
print('\n(3) coherent M1 radiation reaction on the transverse components (reduced equations):')
def run(g1,g2,dpsi0,steps=200000,eta=1e-3):
    # complex transverse spin amplitudes z_i = t_i e^{i psi_i}; rotation at omega_f = 1; moments mu_i = g_i z_i;
    # radiation reaction field B_rr ~ -eta * d/dt(M) with M = sum g_i z_i rotating: d^3M/dt^3 = -i^3 ... -> B_rr = eta * i * M (in the rotating frame, i.e. 90 deg ahead)
    # torque on spin i: dS_i/dt = g_i S_i x B_rr; in the transverse plane its dissipative part reduces z_i along M: dz_i/dt = -eta g_i M (secular, rotating frame)
    z1=complex(t,0); z2=t*complex(math.cos(dpsi0),math.sin(dpsi0)); h=0.05
    for _ in range(steps):
        M=g1*z1+g2*z2; z1+=-eta*g1*M*h; z2+=-eta*g2*M*h
    return z1,z2
random.seed(481)
for name,g1,g2 in (('Ps (e-: -1, e+: +1)',-1.0,1.0),('H (e-: -1, p: +0.0015)',-1.0,0.00152)):
    out=[]
    for k in range(6):
        d0=random.uniform(0,2*math.pi); z1,z2=run(g1,g2,d0)
        dpsi=(math.atan2(z2.imag,z2.real)-math.atan2(z1.imag,z1.real))%(2*math.pi)
        out.append(f'start {d0:.2f} -> dpsi {dpsi:.3f}, t1 {abs(z1):.3f} t2 {abs(z2):.3f}')
    print(f'   {name}:'); [print('      '+o) for o in out]
print('   (Ps: the common mode g1 z1 + g2 z2 decays; spins end in phase with equal reduced transverse amplitude -- one locked phase, 0;\n    H: the electron transverse component decays, the proton keeps its phase -- no locking)')
