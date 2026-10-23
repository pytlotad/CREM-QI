# Audit 592: angular momentum and energy carried by each engine photon (CREM_EMISSION_LEDGER lines).
import glob, re
hbar=1.054571817e-34
rows=[]
for f in sorted(glob.glob('../section-585/runs/*.err')+glob.glob('../section-589/runs/*.err')):
    for line in open(f):
        if not line.startswith('CREM_LEDGER t'): continue
        t=line.split(); eg=float(t[4]); harm=t[6]
        b=t.index('B'); a=t.index('A')
        Eb,Lb=float(t[b+2]),float(t[b+4]); Ea,La=float(t[a+2]),float(t[a+4])
        rows.append((f.split('/')[-1],harm,Lb/hbar,La/hbar,(La-Lb)/hbar,(Eb-Ea-eg)/eg))
worst=0
print("file                 harm  L_before  L_after   Delta L/hbar   (E_b - E_a - Eg)/Eg")
for r in rows:
    print(f"{r[0]:20s} {r[1]:>4s}  {r[2]:8.5f}  {r[3]:8.5f}  {r[4]:+10.6f}   {r[5]:+.2e}")
    worst=max(worst,abs(abs(r[4])-1))
print(f"photons {len(rows)};  max | |Delta L|/hbar - 1 | = {worst:.2e}")
