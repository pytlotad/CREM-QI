# Audit 512: photon as a fast event (current configuration) -- A: invariants before/after each photon; B: cost of snapping |L'| (and |J'|)
# to the Langer value with a tilted photon of |J_gamma| = hbar.  python3 ana512.py <dir>
import math, sys, glob, os
hbar=1.054571817e-34; eV=1.602176634e-19
def V(a): return [float(x) for x in a[:3]]
def dot(a,b): return sum(x*y for x,y in zip(a,b))
def nrm(a): return math.sqrt(dot(a,a))
def unit(a): n=nrm(a); return [x/n for x in a]
def add(a,b): return [x+y for x,y in zip(a,b)]
def mul(a,k): return [x*k for x in a]
def frac(x): return x-math.floor(x)
def fd(a,b): d=abs(frac(a)-frac(b)); return min(d,1-d)
def ks_u(x):
    x=sorted(x); n=len(x); d=max(max((i+1)/n-v,v-i/n) for i,v in enumerate(x)); lam=(math.sqrt(n)+0.12+0.11/math.sqrt(n))*d
    return d,max(0,min(1,2*sum((-1)**(k-1)*math.exp(-2*k*k*lam*lam) for k in range(1,100))))
def load(files):
    out=[]; B=None
    for f in files:
        for l in open(f):
            if l.startswith('BFIELD'): B=[float(x) for x in l.split()[3:6]]
            if not l.startswith('CREM_LEDGER '): continue
            t=l.split(); b=[k for k,x in enumerate(t) if x=='|']; d={'Bf':B,'sph':V(t[t.index('sphoton')+1:])}
            for side,st in (('B',b[0]),('A',b[1])):
                s=t[st+2:]; g=lambda k: s[s.index(k)+1:]
                d[side]=dict(E=float(g('E')[0]),L=float(g('L')[0]),ax=V(g('ax')),m1=V(g('m1')),m2=V(g('m2')))
            out.append(d)
    return out
d=sys.argv[1]
for name,pat,Bry,sg in (('H','H.*.err',13.5983,(1,-1)),('Ps','Ps.*.err',6.80285,(-1,1)),('H 1 T','H1T.*.err',13.5983,(1,-1))):
    ph=load(sorted(glob.glob(os.path.join(d,pat)))); n=len(ph)
    if not n: print(name,'no data'); continue
    def spin(s): return add(mul(unit(s['m1']),sg[0]*0.5),mul(unit(s['m2']),sg[1]*0.5))
    w=dict(dL=0,fL=0,fn=0,S=0,SL=0,fB=0); fJ=[]; fLa=[]; refl=0
    costs=[]; costsJ=[]
    for p in ph:
        b,a=p['B'],p['A']; Lb=b['L']/hbar; La=a['L']/hbar
        if La>Lb: refl+=1
        Lbv=mul(unit(b['ax']),Lb); Lav=mul(unit(a['ax']),La)
        w['dL']=max(w['dL'],abs(nrm(add(Lav,mul(Lbv,-1)))-1)); w['fL']=max(w['fL'],min(fd(Lb,La),fd(Lb,-La)))
        nb=math.sqrt(Bry*eV/abs(b['E'])); na=math.sqrt(Bry*eV/abs(a['E'])); w['fn']=max(w['fn'],fd(nb,na))
        Sb,Sa=spin(b),spin(a); w['S']=max(w['S'],abs(nrm(Sb)-nrm(Sa))); w['SL']=max(w['SL'],abs(dot(Sb,unit(b['ax']))-dot(Sa,unit(a['ax']))))
        if p['Bf'] and nrm(p['Bf'])>0:
            bh=unit(p['Bf']); w['fB']=max(w['fB'],fd(dot(Lbv,bh),dot(Lav,bh)))
        Ja=nrm(add(Lav,Sa)); fJ.append(frac(Ja-0.5)); fLa.append(frac(La-0.5))
        # B: snap |L'| to the Langer value with a tilted hbar photon
        lo=abs(Lb-1); tgt=math.floor(lo-0.5)+0.5   # |L'| reachable range [|Lb - 1|, Lb + 1] (reflection branch when Lb < 1)
        if tgt<lo-1e-12: tgt+=1
        tgt=max(tgt,0.5)
        if tgt<=Lb+1:
            c=(Lb*Lb+1-tgt*tgt)/(2*Lb); costs.append((math.degrees(math.acos(max(-1,min(1,c)))),tgt-Lb,tgt>Lb,c))
        Jb=nrm(add(Lbv,Sb)); loJ=abs(Jb-1); tJ=math.floor(loJ-0.5)+0.5
        if tJ<loJ-1e-12: tJ+=1
        tJ=max(tJ,0.5)
        cJ=(Jb*Jb+1-tJ*tJ)/(2*Jb); costsJ.append(math.degrees(math.acos(max(-1,min(1,cJ)))))
    print(f'== {name}: photons {n} (reflection L increase {refl})')
    print(f'   A: max ||Delta L|/hbar - 1| {w["dL"]:.2e}; max change frac(L/hbar) {w["fL"]:.2e}; frac(n_eff) {w["fn"]:.2e}; |Delta|S|| {w["S"]:.2e}; |Delta S.L^| {w["SL"]:.2e}'
          +(f'; max change frac(L.B^/hbar) {w["fB"]:.2e}' if 'T' in name else ''))
    dk,pk=ks_u(fLa); dj,pj=ks_u(fJ)
    print(f'   A: after photons frac(|L|/hbar - 1/2) KS vs uniform p {pk:.3f}; frac(|J|/hbar - 1/2) KS p {pj:.3f} (uniform = not set)')
    th=[c[0] for c in costs]; cs=[c[3] for c in costs]
    print(f'   B (snap |L\'| to Langer): tilt theta mean {sum(th)/len(th):.1f} deg, max {max(th):.1f}; Delta|L| in [{min(c[1] for c in costs):+.2f}, {max(c[1] for c in costs):+.2f}] hbar (photon rule -1); L increases in {sum(c[2] for c in costs)}/{len(costs)}; cos(J_g, L^) min {min(cs):+.3f}')
    # compare required cos with the E1 directional-spin law density (3/8)(1+c)^2 on [-1,1]: CDF F(c) = ((1+c)^3)/8
    xs=sorted(cs); nn=len(xs); D=max(max((i+1)/nn-((1+x)**3)/8,((1+x)**3)/8-i/nn) for i,x in enumerate(xs))
    print(f'   B: required cos(J_gamma, L^) vs E1 directional-spin law (3/8)(1+c)^2: KS D {D:.3f} (n {nn}); share with c < 0 required 0 vs law {1/8:.3f}')
    thJ=costsJ; print(f'   B (snap |J\'|): tilt theta mean {sum(thJ)/len(thJ):.1f} deg, max {max(thJ):.1f}')
