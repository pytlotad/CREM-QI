import re, sys, math
files = sys.argv[1:]
data = {}
for f in files:
    case = "E1"
    for line in open(f):
        m = re.match(r"\s{2}(\S.*?) \(os", line)
        if line.startswith("(i)"): case = "E1"
        if m: case = m.group(1)
        m = re.match(r"\s+(.+?)\s+([-+]\d\.\d+) \+- (\d\.\d+)\s+\(oczekiwane ([-+]\d\.\d+), odchylka ([-+]?\d+\.\d+) sigma\)", line)
        if m:
            name = case + " | " + m.group(1).strip()
            data.setdefault(name, []).append((float(m.group(2)), float(m.group(3)), float(m.group(4))))
odd = {k: v for k, v in data.items() if abs(v[0][2]) == 0.0}
chi2 = 0.0; n = 0
print("%-70s %s" % ("wielkosc nieparzysta", "z (ziarna 11 / 22 / 33) | srednia"))
for k, v in odd.items():
    zs = [(a - e) / s for a, s, e in v]
    chi2 += sum(z * z for z in zs); n += len(zs)
    mean = sum(a for a, s, e in v) / len(v); err = math.sqrt(sum(s * s for a, s, e in v)) / len(v)
    print("%-70s %s | %+.2e +- %.1e (%+.2f sigma)" % (k, " ".join("%+.2f" % z for z in zs), mean, err, mean / err))
# chi^2 survival with n dof (Wilson-Hilferty)
x = chi2 / n; zwh = (x ** (1 / 3) - (1 - 2 / (9 * n))) / math.sqrt(2 / (9 * n))
p = 0.5 * math.erfc(zwh / math.sqrt(2))
print("chi^2 = %.1f przy %d stopniach swobody, p = %.2f" % (chi2, n, p))
for k, v in data.items():
    if abs(v[0][2]) > 0:
        print("parzysta %-60s %s" % (k, " ".join("%+.2f" % ((a - e) / s) for a, s, e in v)))
