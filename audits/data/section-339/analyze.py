import math,re,statistics as st,sys,os
D=sys.argv[1]
def kv(l):
    d={}
    for tok in l.split()[1:]:
        if '=' in tok:
            k,v=tok.split('=',1)
            try: d[k]=float(v)
            except: pass
    return d
def traj(prefix,i):
    out=open(f"{D}/{prefix}{i}.out").read().splitlines(); err=open(f"{D}/{prefix}{i}.err").read().splitlines()
    end=[l for l in out if l.startswith('idx')]
    if not end: return None
    tend=float(re.search(r't = ([\d.]+)',end[0]).group(1))*1e-12
    lb=[kv(l) for l in out if l.startswith('CREM_LBAL')]
    ph=[l.split() for l in out if l.startswith('CREM_REACH')]
    e1={}
    for l in err:
        if l.startswith('CREM_E1 '): d=kv(l); e1[d['t']]=d
    # checkpoints before the first photon: LBAL t <= photon checkpoint start
    tph=float(ph[0][14]) if ph else None
    rows=[]
    for d in lb:
        if tph is not None and d['t']>tph: break
        e=e1.get(d['t'])
        if e is None: continue
        rows.append((d['hazardOrbits'],d['Eexpected'],e['force/coulomb'],e['spectral'],e['kinematic/closed'],e['total/coulomb'],e['a']))
    W=sum(r[0] for r in rows)
    avg=lambda j: sum(r[0]*math.log(r[j]) for r in rows)/W
    # photon emission instant ~ start of its checkpoint + hazardOrbits*period of that checkpoint
    last=rows[-1]; 
    return {'t':tend,'lnP':avg(1),'lnF':avg(2),'lnS':avg(3),'lnK':avg(4),'lnT':avg(5),'n':len(rows),'stop':end[0].split('stop ')[1].split(',')[0],
            'photons':int(re.search(r'photons (\d+)',end[0]).group(1))}
res=[]
for i in range(48):
    if not (os.path.exists(f"{D}/p{i}.out") and os.path.exists(f"{D}/o{i}.out")): continue
    p=traj('p',i); o=traj('o',i)
    if p and o: res.append((i,p,o))
print(f"pairs analysed: {len(res)}")
dt=[(o['t']-p['t'])/p['t'] for i,p,o in res]
comp={k:[o[k]-p[k] for i,p,o in res] for k in ('lnP','lnF','lnS','lnK','lnT')}
print(f"measured Delta t/t (ortho - para): mean {st.mean(dt):+.3e}  median {st.median(dt):+.3e}")
print(f"predicted from hazard -Delta<ln P>: mean {-st.mean(comp['lnP']):+.3e}")
for k,name in (('lnF','H2 dipole force in E1 acceleration'),('lnS','H7 spectral factor'),('lnK','H5 kinematic / closed Larmor'),('lnT','total/coulomb (F*S*K check)')):
    print(f"   -Delta<{k}> {name:38s}: mean {-st.mean(comp[k]):+.3e}  sd {st.pstdev(comp[k]):.1e}")
resid=[d+cp for d,cp in zip(dt,comp['lnP'])]
print(f"residual (Delta t/t + Delta<ln P>), i.e. not explained by hazard power: mean {st.mean(resid):+.3e} median {st.median(resid):+.3e}")
from collections import Counter
print("stops para:",Counter(p['stop'] for i,p,o in res)," ortho:",Counter(o['stop'] for i,p,o in res))

# --- period and post-photon segment ---
def periods(prefix,i):
    out=open(f"{D}/{prefix}{i}.out").read().splitlines()
    lb=[kv(l) for l in out if l.startswith('CREM_LBAL')]
    ph=[l.split() for l in out if l.startswith('CREM_REACH')]
    end=[l for l in out if l.startswith('idx')][0]
    tend=float(re.search(r't = ([\d.]+)',end).group(1))*1e-12
    tph=float(ph[0][14])
    Ts=[];W=0;lnT=0
    for a,b in zip(lb,lb[1:]):
        if b['t']>tph: break
        T=(b['t']-a['t'])/a['hazardOrbits']; lnT+=a['hazardOrbits']*math.log(T); W+=a['hazardOrbits']
    # photon instant: last checkpoint before/at tph starts at tph; it covers hazardOrbits orbits of last T
    lastT=math.exp(lnT/W) if W else float('nan')
    hz=[d for d in lb if d['t']==tph]
    tphoton=tph+(hz[0]['hazardOrbits'] if hz else 0)*lastT
    return lnT/W if W else float('nan'), tphoton, tend
dT=[];post=[];pre=[]
for i,p,o in res:
    lp,tpp,tep=periods('p',i); lo,tpo,teo=periods('o',i)
    dT.append(lo-lp); pre.append((tpo-tpp)/tpp); post.append(((teo-tpo)-(tep-tpp))/tep)
print(f"\n-Delta<ln T> (orbital period):            mean {-st.mean(dT):+.3e}  sd {st.pstdev(dT):.1e}")
print(f"photon instant Delta/t (ortho-para):      mean {st.mean(pre):+.3e}")
print(f"post-photon segment contribution /t_end:  mean {st.mean(post):+.3e}")
print(f"hazard prediction incl. period: -(Delta lnP + Delta lnT): mean {-st.mean([a+b for a,b in zip(comp['lnP'],dT)]):+.3e}")

# --- hazard suppression S(e) (table from crem_collapse.hpp eccentricOrbitHazardSuppression) ---
keys=[0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.93,0.95,0.97,0.99]
vals=[1.0,0.995009,0.980140,0.955704,0.922208,0.880335,0.830924,0.774946,0.713472,0.647646,0.578663,0.507742,0.436110,0.364990,0.295604,0.229184,0.167007,0.110485,0.061352,0.036433,0.022231,0.010466,0.000853]
def Sof(e):
    e=min(max(e,0),0.999)
    for i in range(1,len(keys)):
        if e<=keys[i]: t=(e-keys[i-1])/(keys[i]-keys[i-1]); return vals[i-1]+t*(vals[i]-vals[i-1])
    return vals[-1]
def lnS(prefix,i):
    out=open(f"{D}/{prefix}{i}.out").read().splitlines()
    lb=[kv(l) for l in out if l.startswith('CREM_LBAL')]
    ph=[l.split() for l in out if l.startswith('CREM_REACH')]; tph=float(ph[0][14])
    W=0;s=0;e2=0
    for d in lb:
        if d['t']>tph: break
        W+=d['hazardOrbits']; s+=d['hazardOrbits']*math.log(Sof(d['ecc'])); e2+=d['hazardOrbits']*d['ecc']**2
    return s/W, e2/W
dS=[];e2p=[];e2o=[]
for i,p,o in res:
    sp,ep=lnS('p',i); so,eo=lnS('o',i); dS.append(so-sp); e2p.append(ep); e2o.append(eo)
print(f"\n<e^2> weighted: para {st.mean(e2p):.2e}, ortho {st.mean(e2o):.2e}")
print(f"-Delta<ln S(e)> (hazard suppression):      mean {-st.mean(dS):+.3e}  sd {st.pstdev(dS):.1e}")
pred=[-(a+b+c) for a,b,c in zip(comp['lnP'],dT,dS)]
print(f"prediction -(Delta lnP + Delta lnT + Delta lnS): mean {st.mean(pred):+.3e}   measured {st.mean(dt):+.3e}")
r=[x-y for x,y in zip(dt,pred)]
print(f"per-pair residual measured - prediction: mean {st.mean(r):+.3e}, sd {st.pstdev(r):.1e}; correlation(measured,pred) = {st.correlation(dt,pred):.3f}")

# --- same with the exact small-e law f(e) = (1-e^2)^1.5/(1+e^2/2) (S = f to 1e-4 for e <= 0.15, audit 329) ---
f=lambda e:(1-e*e)**1.5/(1+e*e/2)
def lnF(prefix,i):
    out=open(f"{D}/{prefix}{i}.out").read().splitlines()
    lb=[kv(l) for l in out if l.startswith('CREM_LBAL')]
    ph=[l.split() for l in out if l.startswith('CREM_REACH')]; tph=float(ph[0][14])
    W=0;s=0
    for d in lb:
        if d['t']>tph: break
        W+=d['hazardOrbits']; s+=d['hazardOrbits']*math.log(f(d['ecc']))
    return s/W
dF=[lnF('o',i)-lnF('p',i) for i,p,o in res]
print(f"\nwith exact f(e) instead of the linearly interpolated table: -Delta<ln f> = {-st.mean(dF):+.3e}")
print(f"=> predicted ortho-para with correct small-e suppression: {-st.mean([a+b+c for a,b,c in zip(comp['lnP'],dT,dF)]):+.3e}")
mx=max(max(lnS_e for lnS_e in []) if False else 0,0)
