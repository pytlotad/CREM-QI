# Audit 471: spin pictures (vector length, transverse phase correlation) against measurements.
import math, itertools
alpha=7.2973525693e-3; mu0=1.25663706212e-6; hbar=1.054571817e-34; c=299792458.0; muB=9.2740100783e-24
pics={'P0 hbar/2, no cone':(0.5,None),'P1 sqrt(3)/2, cone, free phases':(math.sqrt(0.75),'free'),'P2 sqrt(3)/2, cone, locked phases':(math.sqrt(0.75),'locked'),
      'P3 1 (s+1/2), cone, free phases':(1.0,'free'),'P4 1 (s+1/2), cone, locked phases':(1.0,'locked')}
def S12(l,ph,state):   # <S1.S2>/hbar^2 for state in T+1, T0, S
    if ph is None: return {'T+1':0.25,'T0':0.25,'S':-0.25}[state]          # classical product: parallel / antiparallel vectors
    t2=l*l-0.25                                                              # transverse length^2
    proj={'T+1':0.25,'T0':-0.25,'S':-0.25}[state]
    if ph=='free': trans=0.0
    else: trans={'T+1':0.0,'T0':+t2,'S':-t2}[state]
    return proj+trans
print('O1 hyperfine (21 cm: E(F=1)-E(F=0) = A, ratio 1; triplet degenerate):')
res={}
for name,(l,ph) in pics.items():
    tp,t0,s=S12(l,ph,'T+1'),S12(l,ph,'T0'),S12(l,ph,'S')
    ratio=tp-s; deg=abs(tp-t0)<1e-12
    tot_t=2*l*l+2*tp; tot_s=2*l*l+2*s
    res[name]=dict(O1=(abs(ratio-1)<1e-12 and deg)); 
    print(f'  {name:36s} S1.S2: T+1 {tp:+.3f} T0 {t0:+.3f} S {s:+.3f}; split/A {ratio:.3f}; triplet degenerate {deg}; |S|^2 T {tot_t:.3f} S {tot_s:.3f} (QM 2, 0)')
    res[name]['O5']=(abs(tot_t-2)<1e-12 and abs(tot_s)<1e-12)
# O2 Stern-Gerlach: cone pictures project to +-hbar/2 (moment +-mu_B); P0: continuous projection (457, 470)
for name,(l,ph) in pics.items(): res[name]['O2']=ph is not None
# O3 M1 Bi82+: rotating component of the electron moment about F (vector-model triangle with picture lengths), g = 2
I=4.5; F=5.0; lam=243.87e-9; w=2*math.pi*c/lam; k=math.sqrt(1-(83*alpha)**2); rel=(2*k+1)**2/9; tau_meas=0.351e-3
Aqm=mu0*w**3*muB**2/(3*math.pi*hbar*c**3)*4*I/(2*I+1)*rel
print(f'O3 M1 Bi82+: A_QM {Aqm:.0f} 1/s, measured {1/tau_meas:.0f} 1/s')
for name,(l,ph) in pics.items():
    if ph is None: lS,lI,lF=0.5,I,F     # P0: lengths = projections; upper bound sin = 1
    elif abs(l-1)<1e-9: lS,lI,lF=1.0,I+0.5,F+0.5
    else: lS,lI,lF=math.sqrt(0.75),math.sqrt(I*(I+1)),math.sqrt(F*(F+1))
    if ph is None: sin2=1.0
    else:
        cs=(lF*lF+lS*lS-lI*lI)/(2*lF*lS); sin2=1-cs*cs
    m2=(2*muB*lS)**2*sin2
    R=mu0*w**3*m2/(6*math.pi*hbar*c**3)*rel
    res[name]['O3']=(R/Aqm, 2*R/Aqm, 2*R*tau_meas)
    print(f'  {name:36s} R/A_QM {R/Aqm:.4f}{" (upper bound)" if ph is None else ""}; with x2 (vacuum) {2*R/Aqm:.4f}; x2 vs measured {2*R*tau_meas:.3f}')
# O4 J weights (composite S = 1 length): from audit 461
o4={'P0 hbar/2, no cone':0.253,'P1 sqrt(3)/2, cone, free phases':0.026,'P2 sqrt(3)/2, cone, locked phases':0.026,'P3 1 (s+1/2), cone, free phases':0.0,'P4 1 (s+1/2), cone, locked phases':0.0}
print('O4 J weights, max |deviation| from (2J+1)/(3(2l+1)) (461/462; composite S length: P0 classical sum (production 0.253), P1/P2 sqrt(2), P3/P4 1.5):')
for n in pics: print(f'  {n:36s} {o4[n]:.3f}')
print('\nSUMMARY (yes = matches):')
print(f'  {"picture":36s} O1 21cm  O2 SG  O3 M1(x2 = QM)  O4 J-weights  O5 |S|^2')
for n in pics:
    r=res[n]; o3=r['O3'][1]
    print(f'  {n:36s} {"yes" if r["O1"] else "no ":7s} {"yes" if r["O2"] else "no ":6s} {o3:.3f}{"  yes" if abs(o3-1)<1e-9 else "  no ":6s}    {o4[n]:.3f}{" yes" if o4[n]==0 else (" ~" if o4[n]<0.03 else " no")}      {"yes" if r["O5"] else "no"}')
