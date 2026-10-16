# Audit 483: <S12> in ^3P_J -- full QM, its (m_l, m_S)-diagonal part, and the phase-averaged cone picture.
import math, cmath, random
random.seed(483)
# ---- Clebsch-Gordan via explicit formula (small j) ----
def cg(j1,m1,j2,m2,J,M):
    if m1+m2!=M: return 0.0
    f=math.factorial
    pre=math.sqrt((2*J+1)*f(J+j1-j2)*f(J-j1+j2)*f(j1+j2-J)/f(j1+j2+J+1))
    pre*=math.sqrt(f(J+M)*f(J-M)*f(j1-m1)*f(j1+m1)*f(j2-m2)*f(j2+m2))
    s=0.0
    for k in range(0,20):
        a=[k,j1+j2-J-k,j1-m1-k,j2+m2-k,J-j2+m1+k,J-j1-m2+k]
        if min(a)<0: continue
        s+=(-1)**k/ (f(a[0])*f(a[1])*f(a[2])*f(a[3])*f(a[4])*f(a[5]))
    return pre*s
# ---- l = 1 spherical harmonics and <1m|n_i n_j|1m'> by quadrature ----
def Y1(m,th,ph):
    if m==0: return math.sqrt(3/(4*math.pi))*math.cos(th)
    return (-1 if m==1 else 1)*math.sqrt(3/(8*math.pi))*math.sin(th)*cmath.exp(1j*m*ph)
Nt,Np=120,120; nn={}
for m in (-1,0,1):
    for mp in (-1,0,1):
        acc=[[0j]*3 for _ in range(3)]
        for it in range(Nt):
            x=(it+0.5)/Nt*2-1; th=math.acos(x); w=2/Nt
            for ip in range(Np):
                ph=2*math.pi*(ip+0.5)/Np; wt=w*2*math.pi/Np
                n=(math.sin(th)*math.cos(ph),math.sin(th)*math.sin(ph),math.cos(th))
                f=Y1(m,th,ph).conjugate()*Y1(mp,th,ph)*wt
                for i in range(3):
                    for j in range(3): acc[i][j]+=f*n[i]*n[j]
        nn[(m,mp)]=acc
# ---- two spin-1/2: triplet |S m_S> in the product basis; Pauli matrices ----
sx=[[0,1],[1,0]]; sy=[[0,-1j],[1j,0]]; sz=[[1,0],[0,-1]]; sig=[sx,sy,sz]
up,dn=[1,0],[0,1]
def kron(a,b): return [a[i]*b[j] for i in range(2) for j in range(2)]
trip={1:kron(up,up),0:[(x+y)/math.sqrt(2) for x,y in zip(kron(up,dn),kron(dn,up))],-1:kron(dn,dn)}
def pair_op(i,j):   # sigma1_i sigma2_j on 4-dim
    return [[sig[i][a//2][b//2]*sig[j][a%2][b%2] for b in range(4)] for a in range(4)]
P=[[pair_op(i,j) for j in range(3)] for i in range(3)]
def mel(bra,op,ket): return sum(bra[a].conjugate()*op[a][b]*ket[b] for a in range(4) for b in range(4))
def s1s2(bra,ket): return sum(mel(bra,P[i][i],ket) for i in range(3))
def S12_el(ml,ms,mlp,msp):   # <ml ms| S12 |ml' ms'>
    t=0j
    for i in range(3):
        for j in range(3): t+=3*nn[(ml,mlp)][i][j]*mel([complex(x) for x in trip[ms]],P[i][j],[complex(x) for x in trip[msp]])
    if ml==mlp: t-=s1s2([complex(x) for x in trip[ms]],[complex(x) for x in trip[msp]])
    return t
print('(Q) full QM and its (m_l, m_S)-diagonal part, l = 1, S = 1:')
qfull={}; qdiag={}
for J in (0,1,2):
    M=0 if J==0 else J   # any M; use M = J (or 0 for J = 0)
    comps=[(ml,ms,cg(1,ml,1,ms,J,M)) for ml in (-1,0,1) for ms in (-1,0,1) if ml+ms==M]
    full=sum(a*b*S12_el(ml,ms,mlp,msp) for ml,ms,a in comps for mlp,msp,b in comps).real
    diag=sum(a*a*S12_el(ml,ms,ml,ms) for ml,ms,a in comps).real
    qfull[J]=full; qdiag[J]=diag
    print(f'   J = {J}: <S12> full {full:+.5f}   diagonal part {diag:+.5f}   (import 18 value {[-4,2,-0.4][J]:+.3f})')
# ---- (K) classical cones, n from |Y_1ml|^2, spins per P2, J as CG^2 mixture ----
def sample_n(ml):
    while True:
        x=2*random.random()-1; ph=2*math.pi*random.random(); th=math.acos(x)
        p=abs(Y1(ml,th,ph))**2/(3/(4*math.pi))
        if random.random()<p: return (math.sin(th)*math.cos(ph),math.sin(th)*math.sin(ph),math.cos(th))
def cone_pair(ms):
    t=math.sqrt(2)
    if ms==1: p1,p2,ps1,ps2=1,1,random.uniform(0,2*math.pi),random.uniform(0,2*math.pi)
    elif ms==-1: p1,p2,ps1,ps2=-1,-1,random.uniform(0,2*math.pi),random.uniform(0,2*math.pi)
    else: p1,p2=1,-1; ps1=random.uniform(0,2*math.pi); ps2=ps1          # T0: in phase
    s1=(t*math.cos(ps1),t*math.sin(ps1),p1); s2=(t*math.cos(ps2),t*math.sin(ps2),p2); return s1,s2
def dot(a,b): return sum(x*y for x,y in zip(a,b))
N=200000
print('\n(K) cone picture (P2 phases), J as CG^2 mixture over (m_l, m_S):')
for J in (0,1,2):
    M=0 if J==0 else J
    comps=[(ml,ms,cg(1,ml,1,ms,J,M)**2) for ml in (-1,0,1) for ms in (-1,0,1) if ml+ms==M]
    tot=0.0
    for ml,ms,w in comps:
        acc=0.0
        for _ in range(N//len(comps)):
            n=sample_n(ml); s1,s2=cone_pair(ms); acc+=3*dot(s1,n)*dot(s2,n)-dot(s1,s2)
        tot+=w*acc/(N//len(comps))
    print(f'   J = {J}: <S12> cone {tot:+.4f}   QM diagonal part {qdiag[J]:+.4f}   QM full {qfull[J]:+.4f}   (MC error ~0.01)')
# ---- variant: spin transverse phase locked to the orbital azimuth (n's phi) ----
print('\n(K2) variant -- transverse spin phases locked to the orbital azimuth phi of n (psi = phi + const), J = 0 mixture:')
for const in (0.0,math.pi/2,math.pi):
    J=0; comps=[(ml,ms,cg(1,ml,1,ms,0,0)**2) for ml in (-1,0,1) for ms in (-1,0,1) if ml+ms==0]; tot=0.0
    for ml,ms,w in comps:
        acc=0.0; n_=N//6
        for _ in range(n_):
            n=sample_n(ml); ph=math.atan2(n[1],n[0]); t=math.sqrt(2)
            if ms==0: p1,p2=1,-1
            else: p1=p2=ms
            ps=ph+const; s1=(t*math.cos(ps),t*math.sin(ps),p1); s2=(t*math.cos(ps),t*math.sin(ps),p2)
            acc+=3*dot(s1,n)*dot(s2,n)-dot(s1,s2)
        tot+=w*acc/n_
    print(f'   lock constant {const:.3f}: <S12>(J=0) {tot:+.4f}')
