#!/usr/bin/env python3
# Audit 562: per-photon comparison of two emission ledgers (561 vs 562 binary,
# same seeds) and the booking closures from CREM_DEBUG_HIDDEN_J.
# usage: ana562.py label old.ledger new.ledger new.err
import sys
def photons(path):
    out=[]
    for line in open(path):
        if not line.startswith('CREM_LEDGER t '): continue
        w=line.split()
        out.append((float(w[2]),float(w[4])))   # t, Eg
    return sorted(out)
label,oldp,newp,errp=sys.argv[1:5]
old,new=photons(oldp),photons(newp)
print(f"== {label}: photons old {len(old)} new {len(new)}")
worstT=worstE=0.0; unmatched=0
for t,e in new:
    o=min(old,key=lambda p:abs(p[0]-t))
    dt=abs(o[0]-t)/t; de=abs(o[1]-e)/e
    if dt>1e-2: unmatched+=1; continue
    worstT=max(worstT,dt); worstE=max(worstE,de)
print(f"   matched {len(new)-unmatched}, unmatched (other branch) {unmatched}")
print(f"   worst relative photon time change   {worstT:.3e}")
print(f"   worst relative photon energy change {worstE:.3e}")
book={}; broken=0
for line in open(errp):
    if not line.startswith('HIDDEN_J'): continue
    w=line.split()
    try:
        assert len(w)==12 and w[2]=='t'   # label '|dL+dJhid - dJhid_proc|/hbar' has spaces
        where=w[1]; dj=float(w[5]); clos=float(w[9])
    except (AssertionError,ValueError):
        broken+=1; continue   # line interleaved with another thread's output
    b=book.setdefault(where,[0,0.0,0.0,0.0]); b[0]+=1; b[1]=max(b[1],dj); b[2]+=dj; b[3]=max(b[3],clos)
for where,(n,mx,sm,cl) in sorted(book.items()):
    print(f"   booking {where:9s} n {n:5d}  max |dJhid| {mx:.3e} hbar  sum {sm:.3e} hbar  worst closure {cl:.3e} hbar")
print(f"   interleaved (unparsable) HIDDEN_J lines skipped: {broken}")
if not book: print("   no bookings (Delta J_hid identically zero)")
# C3 for every photon of the NEW ledger, any pair and transition (466's
# ana466.py selects Ps 2 -> 1 by the Ps Rydberg): |Delta(L + S) + J_photon|
# and with the hidden angular momentum, S = sum s_i hbar/2, s_i = sign(q_i) m_i.
import math, importlib.util
spec=importlib.util.spec_from_file_location('a466',__file__.replace('section-562/ana562.py','section-466/ana466.py'))
src=open(spec.origin).read().split('worst={}')[0]
ns={}; exec(src,ns)
w3=w3h=0.0; n=0
for line in open(newp):
    if not line.startswith('CREM_LEDGER t '): continue
    d=ns['parse'](line); B,A=d['B'],d['A']
    Jb=ns['add'](ns['mul'](ns['unit'](B['ax']),B['L']),ns['spinvec'](B))
    Ja=ns['add'](ns['mul'](ns['unit'](A['ax']),A['L']),ns['spinvec'](A))
    r=ns['add'](ns['sub'](Ja,Jb),d['sph'])
    w3=max(w3,ns['norm'](r)/ns['hbar'])
    w3h=max(w3h,ns['norm'](ns['add'](r,ns['sub'](A['hid'],B['hid'])))/ns['hbar']); n+=1
print(f"   C3 all photons (new ledger, n {n}): |dJ + J_ph| {w3:.3e} hbar; with hidden J {w3h:.3e} hbar")
