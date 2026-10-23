# Audit 593: Quigg-Rosner as the l -> 0 limit of the rule-B centrifugal term (atomic-like units: hbar = mu = k = 1, a = 1).
from fractions import Fraction as F
import math
def r2_cl(n,L): return 1.0/(n**3*L)                 # <r^-2> Kepler, a_n = n^2 (a = 1), hbar = 1
def r3_B(n,L):  return 1.0/(n**3*L*(L*L-0.25))      # (mu k)^(3/2)/(a_n^(3/2) L (L^2 - 1/4)) with a_n^(3/2) = n^3
def r3_QM(n,l): return 1.0/(n**3*l*(l+0.5)*(l+1))
worst1=0; worst1q=0
for l in range(1,7):
    for n in range(l+1,13):
        L=l+0.5
        cent=(L*L-0.25)*r3_B(n,L); vp=r2_cl(n,L); qm=l*(l+1)*r3_QM(n,l)
        worst1=max(worst1,abs(cent/vp-1)); worst1q=max(worst1q,abs(cent/qm-1))
print(f"T1  max |(L^2-1/4)<r^-3>_B / <V'>_cl - 1| = {worst1:.1e};  vs QM l(l+1)<r^-3>: {worst1q:.1e}  (l = 1..6, n = l+1..12) -> contact remainder 0")
worst2=0; worst2q=0
for n in range(1,13):
    L=0.5
    lim=min(abs(((L+d)**2-0.25)*r3_B(n,L+d)-r2_cl(n,L))/r2_cl(n,L) for d in (1e-6,1e-8,1e-10))
    nB=r2_cl(n,0.5)/(2*math.pi)                       # mu <V'>/(2 pi hbar^2), k = 1
    nQR=1.0/(2*math.pi*n**3*0.5)                       # engine: hbar/(2 pi a^3 n^3 L)
    psi0=1.0/(math.pi*n**3)                            # QM |psi_ns(0)|^2
    worst2=max(worst2,abs(nB/nQR-1)); worst2q=max(worst2q,abs(nB/psi0-1))
    if n<=3: print(f"T2  n={n}: limit of rule-B centrifugal at L->1/2 vs <V'>: {lim:.1e};  mu<V'>/(2 pi hbar^2) = {nB:.6f}, n_QR = {nQR:.6f}, |psi(0)|^2 = {psi0:.6f}")
print(f"T2  max |n_B/n_QR - 1| = {worst2:.1e}, |n_B/|psi(0)|^2 - 1| = {worst2q:.1e} (n = 1..12)")
print("T3  QM continuation: l(l+1)<r^-3>_QM at l -> 0 =", [f"{(1/(n**3*0.5)):.4f}" for n in (1,2,3)], " = 2 pi |psi(0)|^2 =", [f"{2*math.pi/(math.pi*n**3):.4f}" for n in (1,2,3)])
