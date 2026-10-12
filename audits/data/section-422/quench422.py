# Audit 422: magnetic quenching of Ps -- QM reference (measured inputs) vs model variants (a) classical spin dynamics,
# (b) Born projection on field eigenstates with the model's classical HFS, (c) with twice it.
import math
h=6.62607015e-34; hbar=h/(2*math.pi); muB=9.2740100783e-24; g=2.00231930436
tau_p=125.14e-12; tau_o=142.04e-9; Gp=1/tau_p; Go=1/tau_o
HFS_meas=203.394e9; HFS_cl=58.4e9; HFS_q=2*HFS_cl
def y_of(B,hfs):
    x=2*g*muB*B/(h*hfs); return x/(1+math.sqrt(1+x*x))
def G_T0(B,hfs):
    y=y_of(B,hfs); return (Go+y*y*Gp)/(1+y*y)
def Bhalf(curve):
    lo,hi=1e-4,1e3; target=0.5*(Go+0.5*(Gp+Go))  # half way from Go to the high-field limit (Gp+Go)/2
    target=Go+0.5*((Gp+Go)/2-Go)
    for _ in range(200):
        mid=math.sqrt(lo*hi)
        if curve(mid)<target: lo=mid
        else: hi=mid
    return math.sqrt(lo*hi)
# (a) classical: two unit spins, mutual precession at omega_J (from the classical HFS), Larmor +-Omega about z.
def classical_mean_w(B,hfs=HFS_cl,T_over=400):
    omJ=2*math.pi*hfs            # precession of each spin about the other (classical coupling, period 1/hfs)
    Om=g*muB*B/hbar              # electron Larmor; positron opposite
    s1=[1.0,0.0,0.0]; s2=[1.0,0.0,0.0]
    def cross(a,b): return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
    def f(s1,s2):
        z=[0,0,1]
        d1=[omJ*c for c in cross(s2,s1)]; d2=[omJ*c for c in cross(s1,s2)]
        l1=cross(z,s1); l2=cross(z,s2)
        return [d1[i]+Om*l1[i] for i in range(3)],[d2[i]-Om*l2[i] for i in range(3)]
    wfast=max(omJ,2*Om); dt=0.02/wfast; T=T_over*2*math.pi/min(omJ,max(2*Om,omJ*1e-3)); n=int(T/dt)
    n=min(n,400000); acc=0.0
    for i in range(n):
        a1,a2=f(s1,s2)
        b1,b2=f([s1[k]+0.5*dt*a1[k] for k in range(3)],[s2[k]+0.5*dt*a2[k] for k in range(3)])
        c1,c2=f([s1[k]+0.5*dt*b1[k] for k in range(3)],[s2[k]+0.5*dt*b2[k] for k in range(3)])
        e1,e2=f([s1[k]+dt*c1[k] for k in range(3)],[s2[k]+dt*c2[k] for k in range(3)])
        s1=[s1[k]+dt/6*(a1[k]+2*b1[k]+2*c1[k]+e1[k]) for k in range(3)]; s2=[s2[k]+dt/6*(a2[k]+2*b2[k]+2*c2[k]+e2[k]) for k in range(3)]
        n1=math.sqrt(sum(c*c for c in s1)); n2=math.sqrt(sum(c*c for c in s2)); s1=[c/n1 for c in s1]; s2=[c/n2 for c in s2]
        acc+=0.5*(1-sum(s1[k]*s2[k] for k in range(3)))
    return acc/n
def G_classical(B): w=classical_mean_w(B); return Gp*w+Go*(1-w)
print("B[T]   QM(meas HFS)  (b) Born HFS_cl  (c) Born 2HFS_cl  (a) classical   [tau of T0-like state, ns]")
for B in [0.001,0.01,0.03,0.1,0.3,0.6,1.0,2.0,4.0,8.0]:
    print(f"{B:6.3f} {1e9/G_T0(B,HFS_meas):12.3f} {1e9/G_T0(B,HFS_cl):14.3f} {1e9/G_T0(B,HFS_q):15.3f} {1e9/G_classical(B):14.3f}")
bq=Bhalf(lambda B:G_T0(B,HFS_meas)); bb=Bhalf(lambda B:G_T0(B,HFS_cl)); bc=Bhalf(lambda B:G_T0(B,HFS_q))
print(f"B_1/2: QM {bq:.4f} T; (b) {bb:.4f} T (ratio {bq/bb:.3f}); (c) {bc:.4f} T (ratio {bq/bc:.3f})")
# classical B_1/2 by bisection on the classical curve (coarse)
lo,hi=1e-3,50.0; target=Go+0.5*((Gp+Go)/2-Go)
for _ in range(14):
    mid=math.sqrt(lo*hi)
    if G_classical(mid)<target: lo=mid
    else: hi=mid
print(f"(a) classical B_1/2 ~ {math.sqrt(lo*hi):.3f} T (ratio QM/classical {bq/math.sqrt(lo*hi):.3f}); high-field <w> at 50 T: {classical_mean_w(50.0):.4f}")
