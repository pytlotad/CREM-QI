import glob,math,statistics as st,sys
hbar=1.054571817e-34; alpha=7.2973525693e-3; c=299792458.0; me=9.1093837015e-31
k=alpha*hbar*c; mu=me/2
def omega_of(L_hbar,e):
    L=L_hbar*hbar; E=-(mu*k*k/(2*L*L))*(1-e*e)
    a=k/(2*abs(E)); return math.sqrt(k/(mu*a**3))
rows=[]; checks=[]
for f in sorted(glob.glob(sys.argv[1]+'/traj-*.txt'),key=lambda s:int(s.split('-')[-1][:-4])):
    idx=int(f.split('-')[-1][:-4]); Lcl=Ecl=Lph=Eph=0.0; nph=0; es=[]; last=None; per_photon=[]
    for line in open(f):
        if line.startswith('CREM_LBAL'):
            d=dict(x.split('=') for x in line.split()[1:]); e=float(d['ecc']); L=float(d['L']); n=int(d['hazardOrbits'])
            Ef=float(d['Eflux']); Lf=float(d['Lflux']); Ex=float(d['Eexpected'])
            w=omega_of(L,e); f_e=(1-e*e)**1.5/(1+e*e/2)
            if Ef>0 and Lf>0: checks.append((e,(Lf*hbar*w)/(Ef*1.0)/(1) ,f_e))
            Lcl+=n*Lf; Ecl+=n*Ex; es.append(e); last=(e,Lf/Ef if Ef>0 else 0)
        elif line.startswith('CREM_REACH'):
            p=[float(x) for x in line.split()[1:]]
            if len(p)<14: continue
            Lb=p[12]; Lr=p[1]; dL=Lb*(1-Lr); Eg=p[4]
            Lph+=dL; Eph+=Eg; nph+=1
            if last and last[1]>0: per_photon.append((last[0],(dL/Eg)/last[1]))
    status=open(f).read().strip().splitlines()[-1]
    if nph>0 and Lcl>0:
        rows.append((idx,st.mean(es),nph,(Lph/Eph)/(Lcl/Ecl),Lph,Lcl*(Eph/Ecl),status.split(',')[1].strip()))
print("R329a classical flux check: Lflux*omega/Eflux vs f(e) (per checkpoint)")
bins={}
for e,r,fe in checks:
    b=min(int(e*10),9); bins.setdefault(b,[]).append(r/fe)
for b in sorted(bins): v=bins[b]; print(f"   e in [{b/10:.1f},{(b+1)/10:.1f}): n={len(v):4d}  median ratio/f(e) = {st.median(v):.4f}  range {min(v):.3f}-{max(v):.3f}")
print("\nR329c per trajectory: (L_ph/E_ph)/(L_cl/E_cl)")
for r in rows: print(f"   idx {r[0]:2d}  <e>={r[1]:.3f}  photons {r[2]}  ratio {r[3]:.4f}   {r[6]}")
if rows:
    v=[r[3] for r in rows]; print(f"   N={len(v)}  median {st.median(v):.4f}  mean {st.mean(v):.4f}")
    lo=[r[3] for r in rows if r[1]<0.3]; hi=[r[3] for r in rows if r[1]>0.8]
    print(f"   <e> < 0.3: {len(lo)} median {st.median(lo) if lo else float('nan'):.4f};  <e> > 0.8: {len(hi)} median {st.median(hi) if hi else float('nan'):.4f}")
