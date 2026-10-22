# Audit 553: exact nonrelativistic dipole recombination cross section into level n (sum over l) of a hydrogen-like pair, in units of its own Bohr radius and binding
# (atomic units with mu = e = hbar = 1), against Kramers' sigma_n = (32 pi/3 sqrt3) alpha^3 a^2 eta^4/(n (eta^2 + n^2)): the bound-free Gaunt factor g_n(K) = sigma_exact/sigma_Kramers.
# Photoionization (length gauge, energy-normalized continuum R_El = sqrt(2/(pi k)) F_l(-1/k, kr)/r): sigma_PI(nl) = (4 pi^2 alpha omega/3) sum_{l'=l+-1} l_>/(2l+1) |<E l'| r |n l>|^2;
# Milne: sigma_rec(nl) = (g_b/g_+) (omega/(c k))^2 sigma_PI(nl), g_b = 2(2l+1), g_+ = 1, c = 1/alpha (checked: alpha_1s(1e4 K) = 1.58e-13 cm^3/s, Osterbrock & Ferland).  Check: sigma_PI(1s) at threshold = 2^9 pi^2 alpha/(3 e^4) = 0.2252 a^2.
import mpmath as mp
mp.mp.dps=25
alpha=mp.mpf(1)/mp.mpf('137.035999084')
def dipole(n,l,k,lp):
    # Gordon: R_nl = N sum_j a_j r^(l+j) e^(-r/n); F_l'(eta, kr) = C (kr)^(l'+1) e^(-ikr) M(l'+1-i eta, 2l'+2, 2ikr), eta = -1/k;
    # int_0^inf r^s e^(-beta r) M(a,b,cr) dr = Gamma(s+1) beta^-(s+1) 2F1(a, s+1; b; c/beta), beta = 1/n + ik, c = 2ik.
    N=mp.sqrt((mp.mpf(2)/n)**3*mp.factorial(n-l-1)/(2*n*mp.factorial(n+l)))
    eta=-1/k
    C=2**lp*mp.e**(-mp.pi*eta/2)*abs(mp.gamma(lp+1+1j*eta))/mp.gamma(2*lp+2)
    beta=mp.mpf(1)/n+1j*k; cc=2j*k; a=lp+1-1j*eta; b=2*lp+2
    tot=0
    for j in range(n-l):
        aj=(-1)**j*mp.binomial(n+l,n-l-1-j)/mp.factorial(j)*(mp.mpf(2)/n)**(l+j)
        sp=l+lp+j+3
        tot+=aj*mp.gamma(sp+1)/beta**(sp+1)*mp.hyp2f1(a,sp+1,b,cc/beta)
    val=mp.sqrt(2/(mp.pi*k))*N*C*k**(lp+1)*tot
    return mp.re(val) if abs(mp.im(val))<1e-10*abs(val) else abs(val)
def sigma_PI(n,l,k):
    omega=k*k/2+mp.mpf(1)/(2*n*n)
    s=0
    for lp in (l-1,l+1):
        if lp<0: continue
        s+=max(l,lp)/mp.mpf(2*l+1)*dipole(n,l,k,lp)**2
    return 4*mp.pi**2*alpha*omega/3*s
def sigma_rec(n,k):
    omega=k*k/2+mp.mpf(1)/(2*n*n)
    return sum(2*(2*l+1)*(omega*alpha/k)**2*sigma_PI(n,l,k) for l in range(n))
def kramers(n,k):
    eta2=1/(k*k); return (32*mp.pi/(3*mp.sqrt(3)))*alpha**3*eta2**2/(n*(eta2+n*n))
print("check sigma_PI(1s) near threshold (k = 1e-3):", mp.nstr(sigma_PI(1,0,mp.mpf('1e-3')),6), "expected", mp.nstr(2**9*mp.pi**2*alpha/(3*mp.e**4),6))
for KoverB in ('0.0015','0.005','0.015','0.05','0.088','0.15','0.3','0.5','1.0'):
    k=mp.sqrt(2*mp.mpf(KoverB)/2)   # K = k^2/2 in Hartree units of the pair; B = 1/2 -> K/B = k^2
    k=mp.sqrt(mp.mpf(KoverB))
    g=[sigma_rec(n,k)/kramers(n,k) for n in range(1,6)]
    print(f"K/B = {KoverB}: g_n = "+' '.join(mp.nstr(x,4) for x in g))
