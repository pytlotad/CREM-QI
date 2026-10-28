"""Audit 641: conservation at the 2p -> 1s photon from the emission ledger (CREM_LEDGER lines), per system, in two pictures.
Classical (vectors): two objects before (particle 1, particle 2), three after (particle 1, particle 2, photon).  The secular estimator carries the pair as
centre of mass (v_cm) + internal state (orbit E, L vector, moments); particle momenta p_i = m_i v_cm +- mu v_rel sum to M v_cm, so the totals are
  E_tot = sum_i m_i c^2 + E_orbit + (U + F) + kinetic c.m. (relativistic: invariant mass M c^2 = sum m_i c^2 + E_int, energy sqrt((Mc^2)^2 + (Pc)^2)),
  P_tot = p_1 + p_2 (+ p_photon), J_tot = L + S_1 + S_2 (+ s_photon); the photon's orbital part R x p_g cancels the recoil's R x dp (same point R).
Quantum (labels): n, l, S, J, |s_photon|, selection rules E1, and E_g against the level difference of the labels.
usage: python3 ana641.py <run dir>"""
import glob, math, re, sys
c=299792458.0; hbar=1.054571817e-34; eV=1.602176634e-19
me=9.1093837015e-31; mp=1.67262192369e-27
mass={"electron,positron":(me,me,-1,+1),"proton,electron":(mp,me,+1,-1)}
def vec(w,k): i=w.index(k); return [float(x) for x in w[i+1:i+4]]
def add(a,b,s=1): return [x+s*y for x,y in zip(a,b)]
def sc(a,s): return [x*s for x in a]
def norm(a): return math.sqrt(sum(x*x for x in a))
def side(w,m1,m2,q1,q2):
    E=float(w[w.index("E")+1]); L=float(w[w.index("L")+1]); ax=vec(w,"ax")
    M1=vec(w,"m1"); M2=vec(w,"m2"); v=vec(w,"vcm"); UF=float(w[w.index("UF")+1])
    S1=sc(M1,math.copysign(1,q1)*0.5*hbar/norm(M1)); S2=sc(M2,math.copysign(1,q2)*0.5*hbar/norm(M2))
    return dict(E=E,UF=UF,Lv=sc(ax,L/norm(ax)),L=L,S1=S1,S2=S2,S=add(S1,S2),v=v)
pairs={"oPs":"electron,positron","pPs":"electron,positron","oH":"proton,electron","pH":"proton,electron"}
d=sys.argv[1]
for sysname,pair in pairs.items():
    m1,m2,q1,q2=mass[pair]; M=m1+m2; mu=m1*m2/M
    k=8.9875517923e9*eV*eV; Ry1=mu*k*k/(2*hbar*hbar)
    print("=== %s (%s)"%(sysname,pair))
    worst=dict(E=0,Erel=0,P=0,J=0); rules=[]; lines=[]
    for f in sorted(glob.glob("%s/%s_*.txt"%(d,sysname))):
        txt=open(f).read(); led=[l for l in txt.splitlines() if l.startswith("CREM_LEDGER t")]
        born=[l for l in txt.splitlines() if l.startswith("BORN")]
        for l in led:
            parts=l.split("|"); w0=parts[0].split(); B=side(parts[1].split(),m1,m2,q1,q2); A=side(parts[2].split(),m1,m2,q1,q2)
            Eg=float(w0[w0.index("Eg")+1]); n_g=vec(w0,"dir"); s_g=vec(w0,"sphoton")
            # classical, non-relativistic energy with c.m. kinetic terms
            kB=0.5*M*norm(B["v"])**2; kA=0.5*M*norm(A["v"])**2
            dE=(B["E"]+B["UF"]+kB)-(A["E"]+A["UF"]+kA)-Eg
            # relativistic four-momentum: invariant masses
            MB=M*c*c+B["E"]+B["UF"]; MA=M*c*c+A["E"]+A["UF"]
            pB=sc(B["v"],MB/(c*c)); pA=sc(A["v"],MA/(c*c)); pg=sc(n_g,Eg/c)
            EB=math.sqrt(MB**2+(norm(pB)*c)**2); EA=math.sqrt(MA**2+(norm(pA)*c)**2)
            dErel=EB-(EA+Eg); dP=norm(add(add(pB,pA,-1),pg,-1))
            JB=add(B["Lv"],B["S"]); JA=add(A["Lv"],A["S"]); dJ=norm(add(add(JB,JA,-1),s_g,-1))
            worst["E"]=max(worst["E"],abs(dE)/Eg); worst["Erel"]=max(worst["Erel"],abs(dErel)/Eg)
            worst["P"]=max(worst["P"],dP/(Eg/c)); worst["J"]=max(worst["J"],dJ/hbar)
            # quantum labels
            nB=math.sqrt(Ry1/abs(B["E"])); nA=math.sqrt(Ry1/abs(A["E"]))
            lB=B["L"]/hbar-0.5; lA=A["L"]/hbar-0.5; SB=norm(B["S"])/hbar; SA=norm(A["S"])/hbar
            jB=norm(JB)/hbar-0.5; jA=norm(JA)/hbar-0.5
            Jlab=re.findall(r"J=(\d)",born[0]) if born else []
            rules.append((round(nB,9),round(nA,9),round(lB,9),round(lA,9),round(SB,6),round(SA,6),round(norm(s_g)/hbar,12),Jlab[0] if Jlab else "-",round(jB,4),round(jA,4)))
            lines.append((Jlab[0] if Jlab else "-",Eg/eV,dE,dErel,dP/(Eg/c),dJ/hbar,kA/eV))
    print(" classical: worst |dE|/E_g (nonrel, incl. c.m. kinetic) %.2e, relativistic %.2e, |dP|/(E_g/c) %.2e, |dJ|/hbar %.2e"%(worst["E"],worst["Erel"],worst["P"],worst["J"]))
    for r,ln in zip(rules,lines):
        print("  J(Born) %s  n %.9f -> %.9f  l %.9f -> %.9f  S %.6f -> %.6f  |s_g|/hbar %.12f  |J|-1/2: %.4f -> %.4f  E_g %.10f eV  dE %+.2e J  dErel %+.2e J  dP %.1e  dJ %.1e  recoil KE %.3e eV"%(r[7],r[0],r[1],r[2],r[3],r[4],r[5],r[6],r[8],r[9],ln[1],ln[2],ln[3],ln[4],ln[5],ln[6]))
