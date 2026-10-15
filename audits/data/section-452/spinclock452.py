# Audit 452: zitterbewegung (spin) clock -- memory of its resonance with Compton-frequency zero-point modes vs the orbital period.
import math
c=299792458.0; hbar=1.054571817e-34; me=9.1093837015e-31; alpha=7.2973525693e-3; re=2.8179403262e-15
wC=me*c*c/hbar                      # Compton angular frequency (zitter: 2 wC)
tau_e=2*re/(3*c)                   # radiation-reaction time 2 e^2/(3 m c^3)
for wz,name in ((wC,"omega_C"),(2*wC,"2 omega_C (zitter)")):
    G=tau_e*wz*wz                   # radiative linewidth of a charged oscillator at wz
    Q=wz/G; tcoh=1/G
    print(f"{name}: {wz:.4e} rad/s, tau_e {tau_e:.4e} s, linewidth {G:.4e} /s, Q {Q:.2f} (3/(2 alpha) = {1.5/alpha:.2f} for omega_C), coherence {tcoh:.4e} s")
    for sysname,mu_over_m in (("hydrogen",1.0/(1+me/1.67262192369e-27)),("positronium",0.5)):
        for n in (1,2,3):
            wo=mu_over_m*alpha**2*me*c*c/(hbar*n**3)   # Kepler orbital angular frequency = 2|E_n|/(hbar n) ... = mu c^2 alpha^2/(hbar n^3)
            T=2*math.pi/wo
            print(f"   {sysname} n={n}: orbital period {T:.4e} s, coherence/period {tcoh/T:.3e}")
# Ps: relative beat of the two clocks -- equal speeds in the CM frame at every instant (equal masses) -> zero for any orbit.
print("Ps e-/e+ clock beat: v_e = v_p at every instant in the CM frame -> d(tau_e - tau_p)/dt = 0 identically")
