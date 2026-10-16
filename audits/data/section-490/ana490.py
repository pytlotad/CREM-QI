# Audit 490: does the photon set any adiabatic invariant to a quantum value?  Ledger before/after each photon.
import math, sys
sys.path.insert(0,'../section-466')
hbar=1.054571817e-34; eV=1.602176634e-19
def V(a): return [float(x) for x in a[:3]]
def dot(a,b): return sum(x*y for x,y in zip(a,b))
def norm(a): return math.sqrt(dot(a,a))
def unit(a): n=norm(a); return [x/n for x in a]
def add(a,b): return [x+y for x,y in zip(a,b)]
def mul(a,k): return [x*k for x in a]
def parse(line):
    t=line.split(); bars=[k for k,x in enumerate(t) if x=='|']; d={}
    for side,start in (('B',bars[0]),('A',bars[1])):
        seg=t[start+2:]; g=lambda key: seg.index(key)
        d[side]=dict(E=float(seg[g('E')+1]),L=float(seg[g('L')+1]),ax=V(seg[g('ax')+1:]),m1=V(seg[g('m1')+1:]),m2=V(seg[g('m2')+1:]))
    return d
def ks_uniform(x):
    x=sorted(x); n=len(x); d=max(max((i+1)/n-v, v-i/n) for i,v in enumerate(x)); lam=(math.sqrt(n)+0.12+0.11/math.sqrt(n))*d
    return d,max(0.0,min(1.0,2*sum((-1)**(k-1)*math.exp(-2*k*k*lam*lam) for k in range(1,100))))
def frac(x): return x-math.floor(x)
for name,f,B,sg in (('Ps','Ps.err',6.80285,(-1,1)),('H','H.err',13.5983,(1,-1))):
    ph=[parse(l) for l in open(f) if l.startswith('CREM_LEDGER ')]
    worst={'dLnorm':0,'fracL':0,'fracn':0,'S':0,'SL':0,'c12':0}; fJ=[]; c12=[]; refl=0; n=len(ph)
    for d in ph:
        b,a=d['B'],d['A']
        Lb=b['L']/hbar; La=a['L']/hbar
        if La>Lb: refl+=1
        dL=norm(add(mul(unit(a['ax']),a['L']),mul(unit(b['ax']),-b['L'])))/hbar
        worst['dLnorm']=max(worst['dLnorm'],abs(dL-1))
        fb,fa=frac(Lb),frac(La); worst['fracL']=max(worst['fracL'],min(abs(fb-fa),abs(fb-(1-fa)),abs(abs(fb-fa)-1)))
        nb=math.sqrt(B*eV/abs(b['E'])); na=math.sqrt(B*eV/abs(a['E'])); worst['fracn']=max(worst['fracn'],min(abs(frac(nb)-frac(na)),1-abs(frac(nb)-frac(na))))
        def spins(s): return add(mul(unit(s['m1']),sg[0]*0.5),mul(unit(s['m2']),sg[1]*0.5))
        Sb,Sa=spins(b),spins(a); worst['S']=max(worst['S'],abs(norm(Sb)-norm(Sa)))
        worst['SL']=max(worst['SL'],abs(dot(Sb,unit(b['ax']))-dot(Sa,unit(a['ax']))))
        cb=dot(unit(b['m1']),unit(b['m2'])); ca=dot(unit(a['m1']),unit(a['m2'])); worst['c12']=max(worst['c12'],abs(cb-ca)); c12.append(ca)
        Ja=norm(add(mul(unit(a['ax']),La),Sa)); fJ.append(frac(Ja-0.5))
    print(f'== {name}: photons {n} (L increased at {refl}: reflection branch)')
    print(f'   H490a  max | |Delta L vector|/hbar - 1 | = {worst["dLnorm"]:.2e};  max change of frac(L/hbar) (allowing reflection) {worst["fracL"]:.2e}')
    print(f'   H490b  max change of frac(n_eff) {worst["fracn"]:.2e}')
    print(f'   H490c  max |Delta|S|| {worst["S"]:.2e} hbar, max |Delta S.L^| {worst["SL"]:.2e} hbar')
    near=sum(1 for v in fJ if v<0.1 or v>0.9)/n; d,p=ks_uniform(fJ)
    print(f'   H490d  frac(|J|/hbar - 1/2) after photons: within 0.1 of 0 {near:.3f} (uniform 0.200 +- {math.sqrt(0.16/n):.3f}); KS p {p:.3f}')
    d,p=ks_uniform([(x+1)/2 for x in c12])
    print(f'   H490e  max |Delta cos(mu_e,mu_p)| {worst["c12"]:.2e}; after photons KS vs uniform p {p:.3f}; within 0.05 of +-1 / +-0.577 / +-0.5: '
          +' / '.join(f'{sum(1 for x in c12 if abs(abs(x)-t)<0.05)/n:.3f}' for t in (1,1/math.sqrt(3),0.5))+' (uniform 0.050 / 0.100 / 0.100)')
