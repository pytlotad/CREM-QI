"""Audit 623, T2: compare the emission ledger and the END record of image runs (P, C, CP) with the identity run.

usage: python3 compare623b.py T2_dir
Ledger fields are printed in decimal (13-16 digits): they are compared as printed.  END fields are hex (%a): bit for bit.
Expected images: C -- moments m1, m2 -> minus, everything else equal; P -- polar (peri, vcm, dir) -> minus, helicity -> minus,
axial (ax, m1, m2, hidJ, sphoton) and scalars equal; CP -- both.  For every field that fails its image, the field is also
tested unchanged ("stays"), which is what a draw made about the axial orbit axis does under P.
"""
import glob
import os
import sys

POLAR = {"peri", "vcm", "dir"}
PSEUDO = {"hel"}
MOMENT = {"m1", "m2"}
VEC3 = {"dir", "sphoton", "ax", "peri", "m1", "m2", "vcm", "hidJ"}


def parse_ledger_line(line):
    w = line.split()[1:]          # drop CREM_LEDGER
    rec, side, j = {}, "", 0
    while j < len(w):
        key = w[j]
        if key == "|":
            side = w[j + 1] + "_"
            j += 2
            continue
        n = 3 if key in VEC3 else 1
        rec[side + key] = [float(x) for x in w[j + 1:j + 1 + n]]
        j += 1 + n
    return rec


def parse(path):
    runs, cur = {}, None
    for line in open(path):
        if line.startswith("START idx"):
            cur = int(line.split()[2])
            runs[cur] = {"ledger": [], "start": None, "end": None}
        elif line.startswith("CREM_LEDGER_START") and cur is not None:
            w = line.split()
            runs[cur]["start"] = {"m1": [float(x) for x in w[2:5]], "m2": [float(x) for x in w[6:9]]}
        elif line.startswith("CREM_LEDGER t") and cur is not None:
            runs[cur]["ledger"].append(parse_ledger_line(line))
        elif line.startswith("END idx"):
            w = line.split()
            i = int(w[2])
            runs[i]["end"] = {"stop": [float(w[4])], "photons": [float(w[6])], "tau": [float.fromhex(w[8])],
                              "born": [float.fromhex(w[10])], "chosen": [float(w[12])],
                              "vcm": [float.fromhex(x) for x in w[14:17]]}
    return runs


def base(key):
    return key.split("_", 1)[1] if "_" in key else key


def image_sign(T, key):
    k = base(key)
    s = 1.0
    if "P" in T and (k in POLAR or k in PSEUDO):
        s = -s
    if "C" in T and k in MOMENT:
        s = -s
    return s


def tally(stats, key, a, b, s):
    st = stats.setdefault(key, [0, 0, 0])
    st[0] += 1
    if all(x == s * y for x, y in zip(a, b)):
        st[1] += 1
    if all(x == y for x, y in zip(a, b)):
        st[2] += 1


def compare(ref, img, T):
    stats, problems = {}, []
    for i, r in ref.items():
        m = img.get(i)
        if m is None or r["end"] is None or m["end"] is None:
            problems.append("run %d incomplete" % i)
            continue
        if len(r["ledger"]) != len(m["ledger"]):
            problems.append("run %d photons %d vs %d" % (i, len(r["ledger"]), len(m["ledger"])))
        for k in ("m1", "m2"):
            tally(stats, "start_" + k, r["start"][k], m["start"][k], image_sign(T, k))
        for la, lb in zip(r["ledger"], m["ledger"]):
            for k, v in la.items():
                tally(stats, k, v, lb[k], image_sign(T, k))
        for k, v in r["end"].items():
            tally(stats, "END_" + k, v, m["end"][k], image_sign(T, k))
    return stats, problems


def main(d):
    for ref_path in sorted(glob.glob(os.path.join(d, "*_id.err"))):
        stem = ref_path[:-len("_id.err")]
        ref = parse(ref_path)
        nph = sum(len(r["ledger"]) for r in ref.values())
        print("%s: %d przebiegow, %d fotonow" % (os.path.basename(stem), len(ref), nph))
        for T in ("P", "C", "CP"):
            path = stem + "_%s.err" % T
            if not os.path.exists(path):
                print("   %-2s brak pliku" % T)
                continue
            stats, problems = compare(ref, parse(path), T)
            ok = [k for k, s in stats.items() if s[1] == s[0]]
            bad = [k for k, s in stats.items() if s[1] < s[0]]
            print("   %-2s pola zgodne z obrazem %d / %d%s" % (T, len(ok), len(stats),
                  ("; problemy: " + "; ".join(problems)) if problems else ""))
            for k in sorted(bad):
                s = stats[k]
                print("      %-12s obraz %d/%d, bez zmiany %d/%d" % (k, s[1], s[0], s[2], s[0]))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "T2")
