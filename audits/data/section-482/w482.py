# Audit 482: J weights from basins with the spin as a sqrt(S(S+1)) cone object (phase-averaged) instead of a rigid 3/2 hbar vector.
import math, random
random.seed(482)
def qm(J,l): return (2*J+1)/(3*(2*l+1))
def nearestJ_JJ1(J2,l):   # J with J(J+1) nearest to J2, clamped to [l-1, l+1]
    best=min(range(max(0,l-1),l+2),key=lambda J:abs(J*(J+1)-J2)); return best
def langerJ(Jn,l): return min(max(int(round(Jn-0.5)),abs(l-1)),l+1)
N=200000; P=64
res={}
for l in (1,2,3):
    Ll=l+0.5; Lq=math.sqrt(l*(l+1))
    cnt={v:{J:0.0 for J in range(l-1,l+2)} for v in ('V1','V2','V3','V4','V5','V6')}
    for _ in range(N):
        ct=2*random.random()-1; st=math.sqrt(1-ct*ct)   # angle between L (z) and the spin axis
        ax=(st,0,ct)
        # V1 rigid
        J=math.sqrt(Ll*Ll+2.25+2*Ll*1.5*ct); cnt['V1'][langerJ(J,l)]+=1
        # V2 / V3: rigid Langer L, spin cone (|S|^2 = 2, projection 1 on the axis), phase average
        J2=Ll*Ll+2+2*Ll*1*ct; cnt['V2'][langerJ(math.sqrt(J2),l)]+=1; cnt['V3'][nearestJ_JJ1(J2,l)]+=1
        # V4: both cones (|L|^2 = l(l+1), projection l on z), independent phases
        J2b=l*(l+1)+2+2*l*1*ct; cnt['V4'][nearestJ_JJ1(J2b,l)]+=1
        # V5: per-phase assignment with the spin cone and rigid Langer L; V6: both cones, per-phase J(J+1)
        e1=(ct,0,-st); e2=(0,1,0); t=1.0   # transverse spin length: sqrt(2 - 1) = 1
        tl=math.sqrt(l*(l+1)-l*l)
        for k in range(P):
            ps=2*math.pi*(k+0.5)/P
            S=[ax[i]*1+t*(math.cos(ps)*e1[i]+math.sin(ps)*e2[i]) for i in range(3)]
            Jv=[S[0],S[1],S[2]+Ll]; cnt['V5'][langerJ(math.sqrt(sum(x*x for x in Jv)),l)]+=1/P
            pl=2*math.pi*((k*7+3)%P+0.5)/P    # an independent orbital phase (decorrelated sampling)
            Lv=[tl*math.cos(pl),tl*math.sin(pl),l]
            Jw=[Lv[i]+S[i] for i in range(3)]; cnt['V6'][nearestJ_JJ1(sum(x*x for x in Jw),l)]+=1/P
    for v in cnt:
        tot=sum(cnt[v].values()); w={J:cnt[v][J]/tot for J in cnt[v]}
        dev=max(abs(w[J]-qm(J,l)) for J in w); res[(l,v)]=(w,dev)
for v in ('V1','V2','V3','V4','V5','V6'):
    line=f'{v}: '
    for l in (1,2,3):
        w,dev=res[(l,v)]; line+=f' l={l}: '+'/'.join(f'{w[J]:.4f}' for J in sorted(w))+f' (max dev {dev:.4f})'
    print(line)
print('QM: '+'  '.join(f'l={l}: '+'/'.join(f'{qm(J,l):.4f}' for J in range(l-1,l+2)) for l in (1,2,3)))
print(f'(MC: N = {N} orientations, {P} phase nodes for V5/V6; statistical error ~{math.sqrt(0.25/N):.4f})')
