# Audit 590: QM expectation values of the classical leading-order Breit-Pauli operators vs the engine's #18 numbers.
import math
# CODATA 2018 (as physical_constants.hpp)
c=299792458.0; hbar=1.054571817e-34; h=2*math.pi*hbar; e=1.602176634e-19; me=9.1093837015e-31; mp=1.67262192369e-27
alpha=7.2973525693e-3; mu0_4pi=1e-7*(1+5.5e-10); ge=2.00231930436256; gp=5.5856946893
k=alpha*hbar*c                     # e^2/(4 pi eps0)
a0=hbar/(me*c*alpha)
def r3(n,l,mu):                    # <r^-3> QM, reduced mass mu
    a=a0*me/mu; return 1.0/(n**3*l*(l+0.5)*(l+1)*a**3)
# --- Ps 2^3P_J: E_J = A <L.S> + B <S12>, A = (2g-1)/2 k hbar^2 <r^-3>/(m^2 c^2), B = (g/2)^2 k hbar^2 <r^-3>/(4 m^2 c^2)
r=r3(2,1,me/2)
A=(2*ge-1)/2*k*hbar**2*r/(me**2*c**2); B=(ge/2)**2*k*hbar**2*r/(4*me**2*c**2)
LS={0:-2,1:-1,2:1}; S12={0:-4,1:2,2:-0.4}
E={J:A*LS[J]+B*S12[J] for J in (0,1,2)}
p10=(E[1]-E[0])/h/1e9; p21=(E[2]-E[1])/h/1e9
# --- H 2P3/2 - 2P1/2: spin-orbit (k/c^2)[(g-1)/(2 m_e^2) + g/(2 m_e m_p)] hbar^2 <r^-3> Delta<L.S>, Delta<L.S> = 3/2
mu=me*mp/(me+mp); rH=r3(2,1,mu)
fs=k/c**2*((ge-1)/(2*me**2)+ge/(2*me*mp))*hbar**2*rH*1.5/h/1e6
# --- H 2P1/2 hyperfine F=1-0: field at the proton B = mu0/4pi (e/(m_e r^3)) [L + (g_e/2)(3(S.r)r - S)], projected on J:
#     <L.J + (g_e/2)(3(S.r)^2 - S.J)>/j(j+1) = (1 + g_e/2)/(3/4); E_F = (g_p mu_N/hbar) B_J.I, Delta<I.J> = hbar^2 (F=1 minus F=0)
muN=e*hbar/(2*mp)
proj=(1+ge/2)/0.75
hfs=gp*muN*mu0_4pi*e/me*rH*proj*hbar/h/1e6    # Delta E = g_p mu_N (mu0/4pi)(e/m_e)<r^-3> proj Delta<I.J>/hbar, Delta<I.J> = hbar^2
eng={'Ps 2^3P1-2^3P0 [GHz]':(5.4852,p10),'Ps 2^3P2-2^3P1 [GHz]':(4.3856,p21),'H 2P3/2-2P1/2 [MHz]':(10968.705,fs),'H 2P1/2 F=1-0 [MHz]':(59.1687,hfs)}
print("quantity                 engine #18     QM <BP>        engine/QM - 1")
for kk,(en,q) in eng.items(): print(f"{kk:24s} {en:12.4f}  {q:12.4f}   {en/q-1:+.2e}")
