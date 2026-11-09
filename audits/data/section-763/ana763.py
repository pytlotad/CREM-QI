# Audit 763: switch vs default -- photon times and energies in p585 hydrogen/Ps probes -> out/ana763.txt
import os, re
D = os.path.dirname(os.path.abspath(__file__)); h = 6.62607015e-34; e = 1.602176634e-19; c = 299792458.0
def led(f):
    return [(float(m.group(1)), float(m.group(2))) for m in (re.search(r'CREM_LEDGER t (\S+) Eg (\S+)', l) for l in open(f)) if m]
ok = True
for t in ['H2p', 'H3p', 'H4p', 'H3d', 'H3s', 'Ps2p']:
    a = led(os.path.join(D, 'out', 'probes_default', t + '.err')); b = led(os.path.join(D, 'out', 'probes_switch', t + '.err'))
    same_t = len(a) == len(b) and all(x[0] == y[0] for x, y in zip(a, b)); dmax = 0.0
    print(f"{t}: fotonow {len(a)}/{len(b)}, czasy {'identyczne' if same_t else 'ROZNE'}")
    for (ta, ea), (tb, eb) in zip(a, b):
        d = (eb - ea) / h / 1e6; dmax = max(dmax, abs(d))
        print(f"   t {ta:.9e} s  E {ea/e:.9f} eV ({h*c/ea*1e9:.6f} nm, {ea/h/1e12:.6f} THz)  przesuniecie {d:+.4f} MHz ({(eb-ea)/e*1e6:+.6f} ueV)")
    if t != 'Ps2p': ok &= same_t and dmax < 2.0
    else: ok &= same_t and dmax == 0.0
print(f"H763c (czasy bez zmian, |dE| < 2 MHz, Ps bez zmian): {'TAK' if ok else 'NIE'}")
