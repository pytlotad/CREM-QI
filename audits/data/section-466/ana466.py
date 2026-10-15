# Audit 466: emission ledger of the Ps 2p -> 1s photon -- classical and quantum conservation laws, per channel.
import sys, math
hbar=1.054571817e-34; c=299792458.0; me=9.1093837015e-31; M=2*me; eV=1.602176634e-19
sgn=(-1.0,1.0)   # electron (first), positron (second): s_i = sign(q_i) m_i/|m_i|
def V(a): return [float(a[0]),float(a[1]),float(a[2])]
def add(a,b): return [a[i]+b[i] for i in range(3)]
def sub(a,b): return [a[i]-b[i] for i in range(3)]
def mul(a,k): return [x*k for x in a]
def dot(a,b): return sum(a[i]*b[i] for i in range(3))
def norm(a): return math.sqrt(dot(a,a))
def unit(a): n=norm(a); return mul(a,1/n)
def parse(line):
    t=line.split(); d={}
    i=t.index('Eg'); d['Eg']=float(t[i+1]); d['harm']=int(t[t.index('harm')+1]); d['hel']=float(t[t.index('hel')+1])
    d['dir']=V(t[t.index('dir')+1:]); d['sph']=V(t[t.index('sphoton')+1:])
    bars=[k for k,x in enumerate(t) if x=='|']
    for side,start in (('B',bars[0]),('A',bars[1])):
        seg=t[start+2:] ; g=lambda key: seg.index(key)
        d[side]=dict(E=float(seg[g('E')+1]),L=float(seg[g('L')+1]),ax=V(seg[g('ax')+1:]),peri=V(seg[g('peri')+1:]),
                     m1=V(seg[g('m1')+1:]),m2=V(seg[g('m2')+1:]),vcm=V(seg[g('vcm')+1:]),UF=float(seg[g('UF')+1]),hid=V(seg[g('hidJ')+1:]))
    return d
def spinvec(s):  # S = (s1 + s2) hbar/2
    return add(mul(unit(s['m1']),sgn[0]*hbar/2),mul(unit(s['m2']),sgn[1]*hbar/2))
def gamma(v): b2=dot(v,v)/c**2; return 1/math.sqrt(1-b2)
def Jlabel(s,l,triplet):
    if not triplet: return l
    S=spinvec(s); tot=add(mul(unit(s['ax']),(l+0.5)*hbar),mul(unit(S),1.5*hbar))
    return min(max(round(norm(tot)/hbar-0.5),abs(l-1)),l+1)
worst={}
def rec(k,v): worst[k]=max(worst.get(k,0.0),abs(v))
for lab in sys.argv[1:]:
    rows=[parse(l) for l in open(lab) if l.startswith('CREM_LEDGER')]
    worst.clear(); counts={}; ex=[]
    for d in rows:
        B,A=d['B'],d['A']; lb=round(B['L']/hbar-0.5); la=round(A['L']/hbar-0.5)
        nb=round(math.sqrt(-6.80285*eV/B['E'])) if B['E']<0 else 0; na=round(math.sqrt(-6.80285*eV/A['E']))
        if not (nb==2 and na==1): counts['other transitions']=counts.get('other transitions',0)+1; continue
        sc=dot(mul(unit(B['m1']),sgn[0]),mul(unit(B['m2']),sgn[1])); triplet=sc>=0
        key='triplet' if triplet else 'singlet'; counts[key]=counts.get(key,0)+1
        Eg=d['Eg']
        # C1 energy (rest frame of the pair before; lab = that frame when vcm_B = 0)
        Wb=(M*c*c+B['E']+B['UF'])*gamma(B['vcm']); Wa=(M*c*c+A['E']+A['UF'])*gamma(A['vcm'])
        rec(key+' C1 energy (Bohr+U+F) residual / Eg',(Wb-Wa-Eg)/Eg)
        WbB=(M*c*c+B['E'])*gamma(B['vcm']); WaB=(M*c*c+A['E'])*gamma(A['vcm'])
        rec(key+' C1b energy Bohr-only residual / Eg',(WbB-WaB-Eg)/Eg)
        rec(key+' Delta(U+F) [eV]',(B['UF']-A['UF'])/eV)
        rec(key+' U+F before [eV]',B['UF']/eV); rec(key+' U+F after [eV]',A['UF']/eV)
        # C2 momentum (kinematics use the Bohr invariant, as the engine does)
        Pb=mul(B['vcm'],gamma(B['vcm'])*(M*c*c+B['E'])/c**2)
        Pa=mul(A['vcm'],gamma(A['vcm'])*(M*c*c+A['E'])/c**2)
        rec(key+' C2 momentum residual / (Eg/c)',norm(add(sub(Pa,Pb),mul(d['dir'],Eg/c)))/(Eg/c))
        # C3 angular momentum
        Jb=add(mul(unit(B['ax']),B['L']),spinvec(B)); Ja=add(mul(unit(A['ax']),A['L']),spinvec(A))
        rec(key+' C3 |dJ + J_photon| [hbar]',norm(add(sub(Ja,Jb),d['sph']))/hbar)
        rec(key+' C3h with hidden J: |dJ + dJ_hid + J_photon| [hbar]',norm(add(add(sub(Ja,Jb),sub(A['hid'],B['hid'])),d['sph']))/hbar)
        rec(key+' |dJ_hidden| [hbar]',norm(sub(A['hid'],B['hid']))/hbar); rec(key+' |J_hidden before| [hbar]',norm(B['hid'])/hbar)
        rec(key+' |J_photon| - hbar [hbar]',(norm(d['sph'])-hbar)/hbar)
        # C4 moments
        rec(key+' C4 |dm1|/|m1|',norm(sub(A['m1'],B['m1']))/norm(B['m1'])); rec(key+' C4 |dm2|/|m2|',norm(sub(A['m2'],B['m2']))/norm(B['m2']))
        # Q1 dl, Q3 dJ, Q4 dM along the axis before, Q5 parity/C
        rec(key+' Q1 dl + 1',la-lb+1)
        Sb=1 if triplet else 0; Jbl=Jlabel(B,lb,triplet); Jal=Sb  # 1s: J = S
        dJ=Jal-Jbl; ok=abs(dJ)<=1 and not (Jbl==0 and Jal==0); counts[key+f' J {Jbl}->{Jal}']=counts.get(key+f' J {Jbl}->{Jal}',0)+1
        rec(key+' Q3 forbidden dJ (count)',0 if ok else 1)
        axb=unit(B['ax']); rec(key+' Q4 dJ.axis + hbar [hbar]',(dot(sub(Ja,Jb),axb)+hbar)/hbar)
        Cb=(-1)**(lb+Sb); Ca=(-1)**(la+Sb); rec(key+' Q5 parity not flipped (count)',0 if (-1)**lb!=(-1)**la else 1)
        rec(key+' Q5 C not conserved (count)',0 if Cb==Ca*(-1) else 1)
        rec(key+' spin |S| [hbar] (classical)',norm(spinvec(B))/hbar)
        rec(key+' axis turn [deg]',math.degrees(math.acos(max(-1,min(1,dot(unit(A['ax']),axb))))))
        rec(key+' L before / hbar - 1.5',B['L']/hbar-1.5); rec(key+' L after / hbar - 0.5',A['L']/hbar-0.5)
        ex.append((Eg/eV,(B['UF']-A['UF'])/eV,Jbl,triplet))
    print(f'== {lab}: photons {len(rows)}; '+', '.join(f'{k} {v}' for k,v in sorted(counts.items())))
    for k in sorted(worst): print(f'   worst {k:55s} {worst[k]:.3e}')
    for t in sorted(set((x[2],x[3]) for x in ex)):
        es=[x for x in ex if (x[2],x[3])==t]
        print(f'   {"triplet" if t[1] else "singlet"} J={t[0]}: n {len(es)}, E_gamma {sum(x[0] for x in es)/len(es):.9f} eV, Delta(U+F) {sum(x[1] for x in es)/len(es)*1e6:.3f} micro-eV')
