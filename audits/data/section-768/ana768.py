# Audit 768: verdicts H768b, H768c, H768d -> out/ana768.txt
import os, re
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out'); h = 6.62607015e-34; e = 1.602176634e-19; c = 299792458.0; ok = True
def sp(f): return [float(m.group(1)) for m in (re.search(r'= ([\d.]+) MHz', x) for x in open(os.path.join(O, f))) if m]
V767 = [10968.715, 3249.988, 1083.328, 1371.089, 457.029, 228.514]; V763d = [10968.704747, 3249.986592, 1083.328864, 1371.088093, 457.029364, 228.514682]
sw, de, osw = sp('p763_switch.txt'), sp('p763_default.txt'), sp('p763_old_switch.txt'); N = ['2P', '3P', '3D', '4P', '4D', '4F']
print("stan: opcja 768 / opcja 763 (stare drzewo) / domyslny [MHz];  768 - 767;  domyslny - 763")
b = True
for i, t in enumerate(N):
    print(f"  {t}: {sw[i]:.6f} / {osw[i]:.6f} / {de[i]:.6f};  {sw[i]-V767[i]:+.6f};  {de[i]-V763d[i]:+.7f}")
    b &= abs(sw[i] - V767[i]) < 0.005 and abs(de[i] - V763d[i]) < 1e-6
print(f"H768b: {'TAK' if b else 'NIE'}")
def led(f): return [(float(m.group(1)), float(m.group(2))) for m in (re.search(r'CREM_LEDGER t (\S+) Eg (\S+)', l) for l in open(f)) if m]
g = True
for t in ['H2p', 'H3p', 'H4p', 'H3s', 'H3d', 'Ps2p']:
    a = led(os.path.join(O, 'probes_default', t + '.err')); s = led(os.path.join(O, 'probes_switch', t + '.err'))
    st = len(a) == len(s) and all(x[0] == y[0] for x, y in zip(a, s)); dm = 0.0
    print(f"{t}: fotonow {len(a)}/{len(s)}, czasy {'identyczne' if st else 'ROZNE'}")
    for (ta, ea), (tb, eb) in zip(a, s):
        d = (eb - ea) / h / 1e6; dm = max(dm, abs(d))
        print(f"   t {ta:.9e} s  E {ea/e:.9f} eV ({h*c/ea*1e9:.6f} nm, {ea/h/1e12:.6f} THz)  przesuniecie {d:+.5f} MHz ({(eb-ea)/e*1e6:+.7f} ueV)")
    g &= st and (dm == 0.0 if t == 'Ps2p' else dm < 0.05)
print(f"H768c: {'TAK' if g else 'NIE'}")
M = 10969.045
print(f"H768d: 2P wobec pomiaru {M}: domyslny {de[0]-M:+.3f}, opcja 763 {osw[0]-M:+.3f}, opcja 768 {sw[0]-M:+.3f} MHz; opcja 768 = {sw[0]*1e6*h/e*1e6:.4f} ueV = {c/sw[0]/1e6*100:.4f} cm")
