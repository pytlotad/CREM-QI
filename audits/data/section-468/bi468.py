# Audit 468 part B: classical M1 spin-flip rate of the model (precessing moments, P/(hbar omega)) vs measurement (209Bi82+)
# and vs the QM formula (hydrogen 21 cm, no measured A).
import math
hbar=1.054571817e-34; c=299792458.0; me=9.1093837015e-31; alpha=7.2973525693e-3; mu0=1.25663706212e-6; muB=9.2740100783e-24
def A_unit(w): return mu0*w**3*muB**2/(3*math.pi*hbar*c**3)        # = alpha w^3 hbar^2/(3 m^2 c^4)
def A_qm(w,I,Z): k=math.sqrt(1-(Z*alpha)**2); return A_unit(w)*4*I/(2*I+1)*(2*k+1)**2/9
# classical, precession about F at the transition frequency, moment component m_perp: R = mu0 w^3 m_perp^2/(6 pi hbar c^3)
def R_cl(w,mperp2,Z): k=math.sqrt(1-(Z*alpha)**2); return mu0*w**3*mperp2/(6*math.pi*hbar*c**3)*(2*k+1)**2/9
# Bi82+
lam=243.87e-9; w=2*math.pi*c/lam; I=4.5; Z=83; tau_meas=0.351e-3; dtau=0.016e-3
Aq=A_qm(w,I,Z); Rmax=R_cl(w,muB**2,Z)                                # m_perp <= |mu_e| = mu_B (g ~ 2 absorbed in the bound factor)
mq2=muB**2*4*I/(2*I+1)                                               # vector model with sqrt(s(s+1)), sqrt(I(I+1)) lengths
Rvm=R_cl(w,mq2,Z)
print(f'Bi82+: omega {w:.4e} rad/s, A_QM {Aq:.1f} 1/s (tau {1e3/Aq:.3f} ms); measured tau {tau_meas*1e3:.3f}({dtau*1e6:.0f}) ms -> A {1/tau_meas:.0f} 1/s')
print(f'  model hbar/2, upper bound (m_perp = mu_e, precession at omega): R <= {Rmax:.1f} 1/s, tau >= {1e3/Rmax:.3f} ms = {tau_meas/(1/Rmax):.2f}^-1 x ... tau_model/tau_meas >= {(1/Rmax)/tau_meas:.2f}; R_max/A_QM = {Rmax/Aq:.4f} (=(2I+1)/(8I) {(2*I+1)/(8*I):.4f})')
print(f'  vector model with QM lengths: R = {Rvm:.1f} 1/s, tau {1e3/Rvm:.3f} ms; R/A_QM = {Rvm/Aq:.4f}; tau/tau_meas = {(1/Rvm)/tau_meas:.2f}')
# hydrogen 21 cm (no measured A; QM formula as reference, labelled)
w21=2*math.pi*1420.405751768e6; Ah=A_qm(w21,0.5,1)
best=max((math.cos(t/2)**4*math.sin(t/2)**2,t) for t in [i*math.pi/10000 for i in range(10001)])
print(f'H 21 cm: A_QM(formula) {Ah:.4e} 1/s; model hbar/2 (two hbar/2 spins at angle th, precession at omega cos(th/2), photon hbar omega):')
print(f'  R(th)/A_QM = cos^4(th/2) sin^2(th/2)/2, max {best[0]/2:.4f} at th = {math.degrees(best[1]):.1f} deg; QM-length vector model: 0.5000')
