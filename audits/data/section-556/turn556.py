# Audit 556: share of the classically forbidden region r < r_t in the quantum radial dipole integral <E l_in | r | n l'> (atomic units of the pair).
import mpmath as mp
mp.mp.dps=20
def Rnl(n,l,r):
    rho=2*r/n; norm=mp.sqrt((mp.mpf(2)/n)**3*mp.factorial(n-l-1)/(2*n*mp.factorial(n+l)))
    return norm*mp.e**(-rho/2)*rho**l*mp.laguerre(n-l-1,2*l+1,rho)
def REl(k,l,r): return mp.sqrt(2/(mp.pi*k))*mp.coulombf(l,-1/k,k*r)/r
def integrand(n,lp,k,lin): return lambda r: REl(k,lin,r)*r*Rnl(n,lp,r)*r*r
def rturn(K,L):   # K + 1/r = L^2/(2 r^2) -> K r^2 + r - L^2/2 = 0
    return (-1+mp.sqrt(1+2*K*L*L))/(2*K)
KoverB=mp.mpf('0.18'); K=KoverB/2; k=mp.sqrt(2*K)
for lin,n,lp in ((0,2,1),(1,1,0),(1,2,0),(1,3,2),(2,2,1),(2,3,1)):
    L=lin+mp.mpf(1)/2; rt=rturn(K,L); f=integrand(n,lp,k,lin); top=40*n*n
    inner=mp.quad(f,[0,rt]); total=mp.quad(f,[0,rt,n*n,4*n*n,12*n*n,top])
    absin=mp.quad(lambda r: abs(f(r)),[0,rt]); abstot=mp.quad(lambda r: abs(f(r)),[0,rt,n*n,4*n*n,12*n*n,top])
    print(f"l_in={lin} -> n={n} l'={lp}: r_t = {mp.nstr(rt,4)} a; integral from r < r_t / total = {mp.nstr(inner/total,4)}; |integrand| share below r_t = {mp.nstr(absin/abstot,4)}")
