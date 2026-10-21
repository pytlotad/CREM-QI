# Audit 544 (as closure541, which was run inline): Ps per-photon closure E_before - E_after - Eg + (UF_before - UF_after) - Eg^2/(2 M c^2), M = 2 m_e
# (no Larmor term: Ps has no orbital moment).  python3 closure544.py files...
import sys
h=6.62607015e-34; me=9.1093837015e-31; c=299792458.0; eV=1.602176634e-19
for f in sys.argv[1:]:
    idx=None
    for l in open(f):
        if l.startswith('START'): idx=int(l.split()[2])
        elif l.startswith('CREM_LEDGER t'):
            t=l.split(); b=[k for k,x in enumerate(t) if x=='|']; sb=t[b[0]+2:]; sa=t[b[1]+2:]; g=lambda s,k: float(s[s.index(k)+1])
            Eg=float(t[t.index('Eg')+1]); rec=Eg*Eg/(2*2*me*c*c)
            cl=(g(sb,'E')-g(sa,'E')-Eg+(g(sb,'UF')-g(sa,'UF'))-rec)/h*1e-6
            print(f"idx {idx}: t {float(t[2])*1e9:.4f} ns, Eg {Eg/eV:.6f} eV, closure after recoil {cl:+.3f} MHz")
