"""Audit 623, T1: compare image trajectories (P, C, CP) with the identity run, bit for bit and relatively.

usage: python3 compare623.py T1_dir
Each file T1/<sys>_B<b>_<T>.txt holds lines 'F idx k t <hex> r1 x y z r2 ... Emech .. Erad .. rad .. E1 .. E2 ..'.
Prad is measured against E_rad/c.  Expected images: C -- everything equal except m1, m2 -> -m1, -m2; P -- polar vectors (r1, r2, P, Prad) -> minus, axial
(m1, m2, J, Jrad) and scalars unchanged; CP -- both.
"""
import glob
import os
import sys

VECS = ["r1", "r2", "m1", "m2", "P", "J", "Prad", "Jrad"]
SCAL = ["t", "Emech", "Erad", "rad", "E1", "E2"]
POLAR = {"r1", "r2", "P", "Prad"}
MOMENT = {"m1", "m2"}


def parse(path):
    runs = {}
    for line in open(path):
        if not line.startswith("F "):
            continue
        w = line.split()
        idx, k = int(w[1]), int(w[2])
        rec, j = {}, 3
        while j < len(w):
            key = w[j]
            if key in VECS:
                rec[key] = [float.fromhex(x) for x in w[j + 1:j + 4]]
                j += 4
            else:
                rec[key] = [float.fromhex(w[j + 1])]
                j += 2
        if line.endswith("\n") and all(key in rec for key in VECS + SCAL):
            runs[(idx, k)] = rec
    return runs


def sign(T, key):
    s = 1.0
    if "P" in T and key in POLAR:
        s = -s
    if "C" in T and key in MOMENT:
        s = -s
    return s


def compare(ref, img, T):
    total = same = 0
    worst = {}
    scale = {}
    perkey = {}
    for rec in ref.values():
        for key, v in rec.items():
            scale[key] = max(scale.get(key, 0.0), max(abs(x) for x in v))
    # Prad nearly cancels (the pair radiates almost no net momentum), so its natural scale is E_rad/c.
    scale["Prad"] = max(scale.get("Prad", 0.0), scale.get("Erad", 0.0) / 299792458.0)
    missing = len(set(ref) ^ set(img))
    for key_rk, rec in ref.items():
        if key_rk not in img:
            continue
        other = img[key_rk]
        for key, v in rec.items():
            s = sign(T, key)
            for a, b in zip(v, other[key]):
                total += 1
                pk = perkey.setdefault(key, [0, 0])
                pk[0] += 1
                if a == s * b:
                    same += 1
                    pk[1] += 1
                d = abs(a - s * b) / scale[key] if scale[key] > 0 else 0.0
                worst[key] = max(worst.get(key, 0.0), d)
    return total, same, worst, missing, perkey


def main(d):
    for ref_path in sorted(glob.glob(os.path.join(d, "*_id.txt"))):
        base = ref_path[:-len("_id.txt")]
        ref = parse(ref_path)
        print("%s: %d ramek" % (os.path.basename(base), len(ref)))
        for T in ("P", "C", "CP"):
            path = base + "_%s.txt" % T
            if not os.path.exists(path):
                print("   %-2s brak pliku" % T)
                continue
            total, same, worst, missing, perkey = compare(ref, parse(path), T)
            top = sorted(worst.items(), key=lambda kv: -kv[1])[:3]
            print("   %-2s bit w bit %d / %d, ramek bez pary %d, max odchylka wzgledna %.2e (%s)" % (
                T, same, total, missing, max(worst.values()) if worst else 0.0,
                ", ".join("%s %.1e" % kv for kv in top)))
            if same < total:
                print("      " + "; ".join("%s %d/%d %.1e" % (k, perkey[k][1], perkey[k][0], worst[k])
                                         for k in VECS + SCAL if k in perkey))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "T1")
