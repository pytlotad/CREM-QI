import sys
me=9.1093837015e-31; mu=me/2; eV=1.602176634e-19
def kv(l): return {k:float(v) for k,v in (x.split('=',1) for x in l.split()[1:])}
def channels(base):
    outs=[];cur={'en':{},'ph':[]}
    for l in open(base+'.out'):
        if l.startswith('idx'): outs.append((l.split(':')[0].split()[2],cur)); cur={'en':{},'ph':[]}
        elif l.startswith('CREM_SPINENERGY'): d=kv(l); cur['en'][d['t']]=d
        elif l.startswith('CREM_REACH'): cur['ph'].append(float(l.split()[14]))
    skips=[];curs={}
    for l in open(base+'.err'):
        if l.startswith('CREM_SKIP'):
            d=kv(l)
            if d['t']==0.0 and curs: skips.append(curs); curs={}
            curs[d['t']]=d
    skips.append(curs)
    return [(n,c,s) for (n,c),s in zip(outs,skips)]
for base in sys.argv[1:]:
    for name,c,sk in channels(base):
        prev=None; worst=0.0; rows=[]
        for t in sorted(c['en']):
            if t not in sk: continue
            en=c['en'][t]; Eb=-sk[t]['Emag']*mu
            dcred=en['credited']-(prev['credited'] if prev else 0.0); Ea=Eb-dcred; tot=Ea+en['U']
            photon=prev is not None and any(prev['t']<=tp<=t for tp in c['ph'])
            if prev and not photon: worst=max(worst,abs(tot-prev['tot'])/abs(tot))
            rows.append((t,Ea/eV,en['U']/eV*1e3,tot/eV,en['pending']/eV*1e3,int(en['clamps']),photon))
            prev={'t':t,'credited':en['credited'],'tot':tot}
        print(f"{base.split('/')[-1]} {name}: worst rel change of E_orb+U between photon-free checkpoints {worst:.2e}; clamps {rows[-1][5]}; final pending {rows[-1][4]:.3e} meV")
        for r in rows: print(f"    t={r[0]*1e12:9.4f} ps{' *photon before*' if r[6] else '               '} E_orb={r[1]:.9f} eV U={r[2]:+.5f} meV E+U={r[3]:.9f} eV pending={r[4]:.2e} meV")
