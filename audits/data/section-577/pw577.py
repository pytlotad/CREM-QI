# Audit 577: midpoint rule (560) for capture partial waves vs exact QM per channel (Gordon integrals).
import math, numpy as np, mpmath as mp
src=open('../section-555/pw555.py').read().split("for KoverB,label in")[0]
exec(src)   # qm_P, cl_P, ns (gaunt553), alpha
def qm_channel(lin,n,lp,kw):
    kw=mp.mpf(kw); omega=kw*kw/2+mp.mpf(1)/(2*n*n)
    D=ns['dipole'](n,lp,kw,lin)
    sPI=4*mp.pi**2*ns['alpha']*omega/3*max(lp,lin)/mp.mpf(2*lp+1)*D**2
    srec=2*(2*lp+1)*(omega*ns['alpha']/kw)**2*sPI
    return float(srec*kw*kw/(mp.pi*(2*lin+1)))
def modes(E,L,omega):
    """c_+ and c_- of the acceleration of a Kepler orbit (energy E, L) at omega; hyperbola (E > 0) or ellipse (E < 0)."""
    if E>0:
        ah=1/(2*E); e=math.sqrt(1+2*E*L*L); T=math.sqrt(ah**3)
        F=np.linspace(-25,25,400001); t=T*(e*np.sinh(F)-F); x=ah*(e-np.cosh(F)); y=ah*math.sqrt(e*e-1)*np.sinh(F); dt=T*(e*np.cosh(F)-1)
        r=np.sqrt(x*x+y*y); ax=-x/r**3; ay=-y/r**3
        cp=np.trapezoid((ax+1j*ay)*np.exp(-1j*omega*t)*dt,F); cm=np.trapezoid((ax-1j*ay)*np.exp(-1j*omega*t)*dt,F)
        return abs(cp)**2, abs(cm)**2, 'continuous'
    # ellipse: periodic; spectrum is a line sum -- equivalent spectral density at omega from the nearest harmonics (report only)
    a=-1/(2*E); e=math.sqrt(max(0.0,1+2*E*L*L)); nmot=a**-1.5; Tp=2*math.pi/nmot
    u=np.linspace(0,2*math.pi,8192,endpoint=False); M=u-e*np.sin(u); t=M/nmot; x=a*(np.cos(u)-e); y=a*math.sqrt(1-e*e)*np.sin(u)
    r=np.sqrt(x*x+y*y); ax=-x/r**3; ay=-y/r**3; dt=(1-e*np.cos(u))/nmot*(2*math.pi/8192)
    s=omega/nmot   # harmonic number (non-integer: spectral density per unit omega over one period ~ |c(omega)|^2 / Tp per period)
    cp=np.sum((ax+1j*ay)*np.exp(-1j*omega*t)*dt); cm=np.sum((ax-1j*ay)*np.exp(-1j*omega*t)*dt)
    return abs(cp)**2, abs(cm)**2, 'period'
def rule_P(lin,n,lp,K,midpoint=False):
    lg=max(lin,lp); omega=K+0.5/n**2
    E=0.5*(K-0.5/n**2) if midpoint else K
    cp2,cm2,kind=modes(E,lg,omega)
    c2=cp2 if lp==lin-1 else cm2
    dEdw=(2*alpha**3/(3*math.pi))*c2/2.0
    w=2*lg/(2*lin+1)
    return w*dEdw/omega*(1.0/n**3)
for KoverB,label in ((0.18,'Ps'),(0.05,'H')):
    K=KoverB/2; kw=math.sqrt(2*K); nmax=4
    print(f"\n{label}: K/B = {KoverB}")
    print("  l_in -> (n, l'):  QM          rule(K)     ratio   rule(midpoint E) ratio")
    S={'QM':0,'R':0,'M':0}; logs={'R':[],'M':[],'549':[],'554':[]}
    for lin in range(0,5):
        Ls=np.linspace(lin+0.01,lin+0.99,25)
        c0=sum(cl_P(L,K,nmax,False)*L for L in Ls)/Ls.sum(); c1=sum(cl_P(L,K,nmax,True)*L for L in Ls)/Ls.sum()
        for n in range(1,nmax+1):
            qsum=0
            for lp in (lin-1,lin+1):
                if lp<0 or lp>n-1: continue
                q=qm_channel(lin,n,lp,kw); r=rule_P(lin,n,lp,K); m=rule_P(lin,n,lp,K,True); qsum+=q
                S['QM']+=(2*lin+1)*q; S['R']+=(2*lin+1)*r; S['M']+=(2*lin+1)*m
                if lin<=3 and q>0: logs['R'].append(abs(math.log(r/q))); logs['M'].append(abs(math.log(m/q)) if m>0 else 10)
                print(f"  {lin} -> ({n},{lp}):  {q:.3e}   {r:.3e}   {r/q if q>0 else float('nan'):6.3f}   {m:.3e}   {m/q if q>0 else float('nan'):6.3f}")
            if lin<=3 and qsum>0:
                logs['549'].append(abs(math.log(c0[n]/qsum)) if c0[n]>0 else 10); logs['554'].append(abs(math.log(c1[n]/qsum)) if c1[n]>0 else 10)
    print(f"  total (2l+1) weighted, l_in <= 4, n <= 4: rule(K)/QM {S['R']/S['QM']:.3f}, rule(midpoint)/QM {S['M']/S['QM']:.3f}")
    print(f"  mean |log(model/QM)| (l_in <= 3; channels for the rule, per (l_in, n) for 549/554): rule(K) {np.mean(logs['R']):.3f}, rule(mid) {np.mean(logs['M']):.3f}, 549 no-hbar {np.mean(logs['549']):.3f}, 554 hbar {np.mean(logs['554']):.3f}")
