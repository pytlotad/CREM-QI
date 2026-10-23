# Audit 596: the E1 channel selection -- classical ellipse existence with rule-Q values vs QM; planar dipole has no m = 0.
import numpy as np
mismatch=0; total=0
for n in range(2,31):
    for l in range(0,n):
        for n2 in range(1,n):
            for l2 in (l-1,l+1):
                if l2<0: continue
                total+=1
                ellipse=(l2+0.5)<=n2          # L' <= L_circ(E') = n' hbar
                qm=l2<=n2-1
                mismatch+=ellipse!=qm
print(f"T1  channels {total}, ellipse-existence vs QM mismatches: {mismatch}")
for e in (0.1,0.66,0.95):
    M=np.linspace(0,2*np.pi,4096,endpoint=False); u=M.copy()
    for _ in range(100): u-=(u-e*np.sin(u)-M)/(1-e*np.cos(u))
    d=np.stack([np.cos(u)-e, np.sqrt(1-e*e)*np.sin(u), np.zeros_like(u)],1)   # orbit in xy, L along z
    # spherical components about z: d_{+-1} = (dx -+ i dy)/sqrt2, d_0 = dz
    d0=np.sqrt(np.mean(d[:,2]**2)); dn=np.sqrt(np.mean(np.sum(d*d,1)))
    print(f"T2  e={e}: |d_0|/|d| = {d0/dn:.1e}")
L=np.array([0,0,0.5]); Lp=L-1.0*L/np.linalg.norm(L)
print(f"T3  l = 0: L = 0,5 hbar z^, vector rule L - hbar L^ = {Lp[2]:+.2f} hbar z^, |L'| = {np.linalg.norm(Lp):.2f} hbar (same l, orbit reversed) -> rule Q: |L'| = 1,5 hbar (l' = 1)")
