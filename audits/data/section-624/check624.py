"""Audit 624: H624b (rule rates printed by the mechanical path = rule functions), H624c (channel histogram of the first
mechanical photon against the rule shares), H624d (legacy switches = binary before the audit), H624e (seeded streams differ).

usage: python3 check624.py C_dir p624r_binary
C_dir/new/*.err|out: runs on the new rules; C_dir/legacy, C_dir/old: the same 10 seeds with the legacy switches / old binary.
"""
import glob
import math
import os
import re
import subprocess
import sys


def gamma_q(a, x):
    """Regularized upper incomplete gamma Q(a, x) (series / continued fraction)."""
    if x <= 0:
        return 1.0
    gln = math.lgamma(a)
    if x < a + 1:
        s, term, ap = 1.0 / a, 1.0 / a, a
        for _ in range(1000):
            ap += 1
            term *= x / ap
            s += term
            if abs(term) < abs(s) * 1e-15:
                break
        return 1.0 - s * math.exp(-x + a * math.log(x) - gln)
    b, c, d = x + 1 - a, 1e300, 1 / (x + 1 - a)
    h = d
    for i in range(1, 1000):
        an = -i * (i - a)
        b += 2
        d = an * d + b
        d = 1e-300 if abs(d) < 1e-300 else d
        c = b + an / c
        c = 1e-300 if abs(c) < 1e-300 else c
        d = 1 / d
        h *= d * c
        if abs(d * c - 1) < 1e-15:
            break
    return math.exp(-x + a * math.log(x) - gln) * h


def runs(err_path):
    out, cur = {}, None
    for line in open(err_path):
        if line.startswith("START idx"):
            cur = int(line.split()[2])
            out[cur] = {"rules": [], "ladder": []}
        elif cur is not None and line.startswith("RULES"):
            out[cur]["rules"].append(dict(kv.split("=", 1) for kv in line.split()[1:] if "=" in kv))
        elif cur is not None and line.startswith("LADDER"):
            out[cur]["ladder"].append(line.rstrip("\n"))
    return out


def expect(out_path):
    e = {}
    for line in open(out_path):
        if line.startswith("EXPECT"):
            w = line.split()
            e[w[1]] = {int(m.group(1)): float(m.group(2)) for m in re.finditer(r"n'=(\d+):(\S+)", line)}
    return e


def chi2(counts, shares):
    keys = sorted(shares)
    total = sum(counts.get(k, 0) for k in keys)
    if total == 0:
        return None
    # merge bins with expected count < 5 into one
    bins, acc_o, acc_e = [], 0, 0.0
    for k in sorted(keys, key=lambda k: shares[k]):
        o, ex = counts.get(k, 0), shares[k] * total
        acc_o += o
        acc_e += ex
        if acc_e >= 5:
            bins.append((acc_o, acc_e))
            acc_o, acc_e = 0, 0.0
    if acc_e > 0:
        if bins:
            o, ex = bins.pop()
            bins.append((o + acc_o, ex + acc_e))
        else:
            bins.append((acc_o, acc_e))
    x = sum((o - ex) ** 2 / ex for o, ex in bins if ex > 0)
    dof = len(bins) - 1
    p = gamma_q(dof / 2.0, x / 2.0) if dof > 0 else float("nan")
    extra = sum(v for k, v in counts.items() if k not in shares)
    return total, x, dof, p, extra


def main(d, rbin):
    states = sorted({re.sub(r"_\d+\.err$", "", os.path.basename(p)) for p in glob.glob(os.path.join(d, "new", "*.err"))})
    printed = []
    for st in states:
        n0, l0 = map(int, re.match(r"n(\d+)l(\d+)", st).groups())
        allruns = {}
        for p in glob.glob(os.path.join(d, "new", st + "_*.err")):
            allruns.update(runs(p))
        ex = expect(glob.glob(os.path.join(d, "new", st + "_*.out"))[0])
        cats = {}
        nophoton = other_state = 0
        first_times = []
        for i, r in sorted(allruns.items()):
            if not r["rules"]:
                nophoton += 1
                continue
            f = r["rules"][0]
            n, L = float(f["n"]), float(f["L"])
            if abs(n - n0) > 1e-3 or abs(L - (l0 + 0.5)) > 1e-3:
                other_state += 1
                continue
            printed.append((st, n, L, f))
            kind = f["type"] if f["type"] != "E1" else ("sE1" if l0 == 0 else ("up" if f["up_photon"] == "1" else "down"))
            cats.setdefault(kind, {}).setdefault(int(f["n_final"]), 0)
            cats[kind][int(f["n_final"])] += 1
            first_times.append(r["ladder"][0] if r["ladder"] else "")
        print("%s: %d przebiegow, bez fotonu %d, pierwszy foton z innego stanu %d" % (st, len(allruns), nophoton, other_state))
        for kind, counts in sorted(cats.items()):
            line = "   %-5s %s" % (kind, " ".join("n'=%d:%d" % kv for kv in sorted(counts.items())))
            if kind in ex:
                res = chi2(counts, ex[kind])
                if res:
                    total, x, dof, p, extra = res
                    line += "  | regula %s | chi^2 %.2f / %d, p = %.3f%s" % (
                        " ".join("%d:%.3f" % (k, v * total) for k, v in sorted(ex[kind].items())), x, dof, p,
                        ", poza regula %d" % extra if extra else "")
            print(line)
        print("   H624e: rozne pierwsze fotony (LADDER) w %d z %d przebiegow" % (len(set(first_times)), len(first_times)))
    # H624b
    inp = "".join("%s %s\n" % (f["n"], f["L"]) for _, _, _, f in printed)
    res = subprocess.run([rbin], input=inp, capture_output=True, text=True, env=dict(os.environ, PAIR="proton,electron")).stdout.split("\n")
    worst = {"sE1": 0.0, "sE2": 0.0, "up": 0.0}
    for (st, n, L, f), line in zip(printed, res):
        w = line.split()
        for key, val in (("sE1", float(w[2])), ("sE2", float(w[3])), ("up", float(w[4]))):
            pv = float(f[key])
            if val > 0 or pv > 0:
                worst[key] = max(worst[key], abs(pv - val) / max(abs(val), 1e-300))
    print("H624b: %d fotonow, max |druk - funkcja|/funkcja: sE1 %.2e, sE2 %.2e, up %.2e" % (len(printed), worst["sE1"], worst["sE2"], worst["up"]))
    # H624d
    same = tot = 0
    for p in glob.glob(os.path.join(d, "legacy", "*.err")):
        q = os.path.join(d, "old", os.path.basename(p))
        a = [l for l in open(p) if not l.startswith("RULES")]
        b = open(q).readlines() if os.path.exists(q) else []
        tot += 1
        same += a == b
    print("H624d: tryb LEGACY (obie zmienne) wobec binarki sprzed audytu: %d / %d plikow stderr identycznych" % (same, tot))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
