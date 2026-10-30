# Audit 651: photon energy from the 2P / 1S field eigenvalues (CREM_PS_EIGEN_PHOTON=1).
import re,sys,math,glob,os
d=sys.argv[1]; eV=1.602176634e-19; h=4.135667696e-15; muB=5.7883818060e-5; g=2.00231930436; Dm=203.3506e9   # Ishida 2014
P={'3P0':5.1017085689,'3P1':5.1017312538,'3P2':5.1017493913,'1P1':5.1017388096}; HFS=5.1025799222-5.1017388096
M={f'{k}->1{"3S" if s==0 else "1S"}':v+s*HFS for k,v in P.items() for s in (0,1)}
def parse(f):
    t=open(f).read(); r={}
    r['E']=float(re.search(r'Eg (\S+)',t).group(1))/eV
    m=re.search(r'eigen photon k=(\d+) w_k=(\S+) w_q=(\S+) jump (\S+) classical (\S+) import (\S+) gapShift (\S+)',t)
    r['eig']=bool(m)
    if m: r.update(k=int(m.group(1)),wk=float(m.group(2)),wq=float(m.group(3)),imp=float(m.group(6))/eV,gs=float(m.group(7))/eV)
    im=re.findall(r'\(1S\) q=\S+ m=(-?\d+) P_S=(\S+) -> (.+)',t)
    r['img']=(im[-1][2]+' m='+im[-1][0]) if im else '-'
    r['born']=[l for l in t.splitlines() if l.startswith('BORN') and 'eigen photon' not in l]
    r['stop']=int(re.search(r'stop (\d+)',t).group(1))
    return r
def runs(B,c,s,N): return [parse(f'{d}/B{B}_{c}_{s}_{i}.err') for i in range(N)]
def BR(B):
    x=2*g*muB*B/h/Dm; return Dm/2*(math.sqrt(1+x*x)-1)*h*1e6, Dm*math.sqrt(1+x*x)*h*1e6, Dm*h*1e6
print('H651a (B = 1e-6 T, seeds 0-15): line vs nearest zero-field multiplet component (650)')
for c in ['O1','P1']:
    R=runs('1e-6',c,'L2',16)
    dev=[min(abs(r['E']-v) for v in M.values())*1e6 for r in R]
    print(f'  {c}: eigen photons {sum(r["eig"] for r in R)}/16, max dev {max(dev):.2e} ueV, |import| max {max(abs(r.get("imp",0)) for r in R)*1e6:.2e} ueV')
for B,N in [('0.1',32),('0.866',16),('1',32)]:
    print(f'=== B = {B} T')
    for c in ['O1','P1']:
        R=runs(B,c,'L2',N); G={}
        for r in R:
            if r['eig']: G.setdefault((round(r['wk']/1e-30),round(r['wq']/1e-30),r['img']),[]).append(r)
        spread=max((max(x['E'] for x in v)-min(x['E'] for x in v) for v in G.values()),default=0)
        imps=[r['imp']*1e6 for r in R if r['eig']]
        print(f'  {c} L2: stops {sorted({r["stop"] for r in R})}, eigen photons {len(imps)}/{N}, groups (w_k, w_q) {len(G)}, max spread in group {spread:.2e} eV;'
              f' import ueV: min {min(imps):.2f} max {max(imps):.2f} mean|.| {sum(map(abs,imps))/len(imps):.2f}')
        # H651d: same w_k, different 1S state
        mix,brm,hfs=BR(float(B)); out=[]
        byk={}
        for (wk,wq,img),v in G.items(): byk.setdefault(wk,{})[img]=v[0]['E']
        for wk,dd in byk.items():
            t0=[e for i,e in dd.items() if i.startswith('triplet image m=0')]; s=[e for i,e in dd.items() if i.startswith("S' image")]
            pm=[e for i,e in dd.items() if i.startswith('triplet image m=1') or i.startswith('triplet image m=-1')]
            for a in pm:
                for b in t0: out.append(('m=+-1 - T0\'',(a-b)*1e6,mix))
            for a in s:
                for b in t0: out.append(("S' - T0'",(a-b)*1e6,brm))
        for name,v,ref in out:
            print(f'      same k: line({name}) = {v:+.4f} ueV; measured-based {ref:.4f} ueV, ratio {v/ref:.5f}')
        if not out: print('      (no k with two different 1S states in this sample)')
for B in ['0.1','1']:
    L=runs(B,'O1','L2',16)
    for s in ['c2','n23']:
        Pp=runs(B,'O1',s,16); same=[i for i in range(16) if Pp[i]['born']==L[i]['born']]
        eq=sum(abs(Pp[i]['E']/L[i]['E']-1)<1e-12 for i in same)
        mx=max([abs(Pp[i]['E']-L[i]['E']) for i in same],default=0)
        print(f'H651c B {B} {s}: labels equal {len(same)}/16, lines equal (1e-12) {eq}/{len(same)}, max |dE| among equal labels {mx:.2e} eV')
print('H651d from w_q (rung, any k): w(T0\') - w(m=+-1) and w(T0\') - w(S\') vs Breit-Rabi from the measured Delta_HFS')
for B,N in [('0.1',32),('0.866',16),('1',32)]:
    W={}
    for c in ['O1','P1']:
        for r in runs(B,c,'L2',N):
            if r['eig']: W.setdefault(r['img'].replace('m=-1','m=1'),set()).add(round(r['wq']/eV*1e6,6))
    mix,brm,_=BR(float(B))
    t0=W.get('triplet image m=0',set()); pm=W.get('triplet image m=1',set()); sp=W.get("S' image m=0",set())
    s=f'  B {B} T: w values (ueV) {dict((k,sorted(v)) for k,v in W.items())}'
    if t0 and pm: s+=f'\n      T0\' - m=+-1 = {min(t0)-max(pm):.4f} ueV, Breit-Rabi {mix:.4f}, ratio {(min(t0)-max(pm))/mix:.5f}'
    if t0 and sp: s+=f'\n      T0\' - S\' = {min(t0)-max(sp):.4f} ueV, Breit-Rabi {brm:.4f}, ratio {(min(t0)-max(sp))/brm:.5f}'
    print(s)
R=runs('1e-6','O1','L2',16)
print(f'1e-6 T o-Ps: seeds with |import| > 1 ueV (k and classical spins give different J): {sum(abs(r["imp"])>1e-6 for r in R)}/16;',
      ' '.join(f'{r["imp"]*1e6:+.2f}' for r in R))
