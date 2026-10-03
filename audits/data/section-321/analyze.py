import sys, math
B=6.80285*1.602176634e-19
for f in sys.argv[1:]:
    rows=[list(map(float,l.split()[1:])) for l in open(f) if l.startswith("CREM_REACH")]
    first=[r for r in rows if r[11]>1.999]
    print(f"{f.split('/')[-1]}: REACH lines {len(rows)}, first photons (n_before > 1.999) {len(first)}")
    if not first: continue
    trimmed=[r for r in first if r[7]<0.9999*r[10]]
    nA=[r[11]/math.sqrt(r[0]) for r in first]; LA=[r[12]*r[1] for r in first]
    print(f"  trimmed {len(trimmed)}/{len(first)}; E_gamma/(3B/4): mean {sum(r[4] for r in first)/len(first)/(0.75*B):.4f}"
          f" min {min(r[4] for r in first)/(0.75*B):.4f} max {max(r[4] for r in first)/(0.75*B):.4f}")
    s=sorted(nA); m=len(s)
    print(f"  n_E after: min {s[0]:.4f} median {s[m//2]:.4f} max {s[-1]:.4f}; at sqrt2 (+-0.01) {sum(abs(x-2**0.5)<0.01 for x in nA)}; below 1.05 {sum(x<1.05 for x in nA)}")
    print(f"  trimmed: max |n_E - L'/hbar| {max([abs(r[11]/math.sqrt(r[0])-r[12]*r[1]) for r in trimmed] or [0]):.2e}, max |e^2 after| {max([abs(r[6]) for r in trimmed] or [0]):.2e}")
    print(f"  untrimmed e^2 after: {[round(r[6],4) for r in first if r not in trimmed]}")
    print(f"  harmonics: {sorted(set(int(r[10]) for r in first))}")
