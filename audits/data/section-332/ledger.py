import sys,math
hbar=1.054571817e-34; me=9.1093837015e-31; mu=me/2; eV=1.602176634e-19
for ph in (1,2):
    f=f"{sys.argv[1]}/ch{ph}.txt"; skips=[]; photons=[]; runs=0
    for l in open(f):
        if l.startswith('CREM_SKIP'):
            d=dict(x.split('=') for x in l.split()[1:]); row={k:float(v) for k,v in d.items()}
            if row['t']==0.0 and skips and not row.get('_'): runs+=1
            if runs==0: skips.append(row)
        elif l.startswith('CREM_REACH'):
            p=[float(x) for x in l.split()[1:]]
            if len(p)>=14: photons.append(p)
    name='para ' if ph==1 else 'ortho'
    # |J| from L (magnitude along Lhat) and S (projection S, magnitude Smag)
    def Jmag(s): Sp=s['S']; Sm=s['Smag']; return math.sqrt((s['Lorb']+Sp)**2+max(Sm*Sm-Sp*Sp,0))
    tph=[p[13] for p in photons]
    segs=[[] for _ in range(len(tph)+1)]
    for s in skips:
        seg=sum(1 for t in tph if s['t']>t); segs[seg].append(s)
    print(f"{name}: checkpoints {len(skips)}, photons {len(photons)} at t = {[round(t*1e12,3) for t in tph]} ps")
    for i,sg in enumerate(segs):
        if len(sg)<2: print(f"   segment {i}: {len(sg)} checkpoints"); continue
        J=[Jmag(s) for s in sg]; E=[s['Emag'] for s in sg]; L=[s['Lorb'] for s in sg]; Sm=[s['Smag'] for s in sg]
        print(f"   segment {i}: {len(sg)} cp, |J| {J[0]:.6f}->{J[-1]:.6f} (max dev {max(abs(x-J[0]) for x in J):.2e}), "
              f"|E| rel change {(E[-1]-E[0])/E[0]:+.2e}, |L| {L[0]:.5f}->{L[-1]:.5f}, |S| {min(Sm):.4f}..{max(Sm):.4f}")
    # whole-trajectory energy ledger (specific energies * mu)
    E0=-skips[0]['Emag']*mu; E1=-skips[-1]['Emag']*mu
    Eg=sum(p[4] for p in photons)
    print(f"   energy: orbital E start {E0/eV:.4f} eV, end {E1/eV:.4f} eV, drop {(E0-E1)/eV:.4f} eV; photons {Eg/eV:.4f} eV; drop - photons = {(E0-E1-Eg)/eV:+.4e} eV (recoil + rest)")
