# Audit 587: the annihilation QED block against photon emission (imports 1, 2), offline.
import numpy as np, math
src=open('../section-560/tempo_e1_punkt_srodkowy.py').read().split('if __name__ == "__main__":')[0]
g={}; exec(src,g)
# constants as in physical_constants.hpp (CODATA 2018)
hbar=1.054571817e-34; c=299792458.0; m=9.1093837015e-31; alpha=7.2973525693e-3
re=alpha*hbar/(m*c)
lhs=2*math.pi*alpha*hbar**3/(m*m*c); rhs=hbar*(4*math.pi*re*re*c)/(2*alpha)
print(f"T2  2 pi alpha hbar^3/(m^2 c) = {lhs:.15e}   hbar sigma v/(2 alpha) = {rhs:.15e}   rel {abs(lhs-rhs)/lhs:.1e}")
cs=np.linspace(-1,1,10001); wpi=0.5*(1-cs*cs)
d=np.max(np.abs((1-cs*cs/3)-(2/3)*(1+wpi)))
print(f"T3  max |(1 - c^2/3) - (2/3)(1 + w_pi)| over 10^4 angles = {d:.1e}")
print("    general even form A|eps|^2 + B|eps.n|^2; m=+-1: |eps.n|^2 = (1-c^2)/2 -> B/A = 1 (ratio normal/in-plane = 2); m=0 -> 1 + c^2;"
      f" m-sum isotropic: {np.max(np.abs(2*(1+0.5*(1-cs*cs))+(1+cs*cs)-4)):.1e}")
# T4: Debye asymptotic of the co-rotating mode amplitude of a Kepler orbit, checked against the 560 Fourier
def logJ(s,e):     # ln J_s(s e) ~ -s(alpha - tanh alpha) - 0.5 ln(2 pi s tanh alpha), cosh alpha = 1/e
    a=math.acosh(1/e); t=math.tanh(a); return -s*(a-t)-0.5*math.log(2*math.pi*s*t)
def logJp(s,e):    # ln J'_s(s e) ~ -s(alpha - tanh alpha) + 0.5 ln(sinh(2 alpha)/(4 pi s))
    a=math.acosh(1/e); t=math.tanh(a); return -s*(a-t)+0.5*math.log(math.sinh(2*a)/(4*math.pi*s))
def mode(s,e):     # |c_s^+|/a = (A + B)/2, A = (2/s) J'_s(se), B = 2 sqrt(1-e^2) J_s(se)/(se)  (432)
    A=(2/s)*math.exp(logJp(s,e)); B=2*math.sqrt(1-e*e)*math.exp(logJ(s,e))/(s*e); return 0.5*(A+B)
e=math.sqrt(3)/2
print("T4  check of the asymptotic against the 560 Fourier (n = 1, L = 1/2): s, |c+|/a Fourier, asymptotic, ratio")
for s in (20,50,100,200):
    cp,_=g['fourier_components'](1.0,0.5,s,N=1<<16)
    print(f"    {s:4d}  {cp:.4e}  {mode(s,e):.4e}  {mode(s,e)/cp:.4f}")
Ry=13.605693122994*0.5            # Ps reduced-mass Rydberg [eV]
hOmega=2*Ry                        # hbar Omega at n = 1
s_ann=511e3*0.99999/hOmega         # one photon of 2 gamma carries ~ m c^2
a=math.acosh(1/e); t=math.tanh(a)
log10A=(-s_ann*(a-t))/math.log(10)
print(f"    Ps 1s: hbar Omega = {hOmega:.3f} eV, s = m c^2/(hbar Omega) = {s_ann:.0f}; mode amplitude ~ 10^{log10A:.0f}, power ~ 10^{2*log10A+4*math.log10(s_ann):.0f} of the fundamental's scale")
