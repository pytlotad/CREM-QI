import glob,statistics as st,math,sys
dl=[];trim=0;nph=0;Lph=Lcl=Eph=Ecl=0.0;byE={}
for f in glob.glob(sys.argv[1]+'/traj-*.txt'):
    lastE=None; tLcl=tEcl=tLph=tEph=0.0; n=0
    for line in open(f):
        if line.startswith('CREM_LBAL'):
            d=dict(x.split('=') for x in line.split()[1:]); h=int(d['hazardOrbits'])
            tLcl+=h*float(d['Lflux']); tEcl+=h*float(d['Eexpected']); lastE=float(d['ecc'])
        elif line.startswith('CREM_REACH'):
            p=[float(x) for x in line.split()[1:]]
            if len(p)<14: continue
            d_=p[12]*(1-p[1]); dl.append(d_); n+=1
            if p[7]<0.9999*p[10]: trim+=1
            e2=p[5]; b=min(int(math.sqrt(max(e2,0))*5),4); byE.setdefault(b,[]).append(d_)
            tLph+=d_; tEph+=p[4]
    if n: Lph+=tLph; Lcl+=tLcl; Eph+=tEph; Ecl+=tEcl; nph+=n
print(f"photons {nph}, trimmed by ceiling {trim} ({100*trim/nph:.0f}%)")
print(f"dL per photon [hbar]: mean {st.mean(dl):.4f} +- {st.stdev(dl)/math.sqrt(len(dl)):.4f}, median {st.median(dl):.4f}, range {min(dl):.3f}..{max(dl):.3f}, negative {sum(x<0 for x in dl)}")
for b in sorted(byE): v=byE[b]; print(f"   e before in [{b/5:.1f},{(b+1)/5:.1f}): n={len(v):3d} mean dL {st.mean(v):.4f} hbar")
print(f"sums over trajectories: L_ph {Lph:.3f} hbar vs L_cl {Lcl:.3f} hbar  -> ratio {Lph/Lcl:.4f};  E_ph/E_cl {Eph/Ecl:.4f}")
