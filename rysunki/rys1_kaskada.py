"""Rys. 1 do Artykul.md: kaskada E1 wodoru (a) i ta sama kaskada na plaszczyznie dzialan (b).

Udzialy QM: wzor Gordona (audits/data/section-560/tempo_e1_punkt_srodkowy.py, rate_qm).
Udzialy silnika (kanaly Delta l = -1): audits/data/section-613/p613.txt.
Stosunki czasow zycia silnik/QM: audyty 560, 612 (Model.md, sekcja 9).
Uruchomienie: python3 rysunki/rys1_kaskada.py  ->  rysunki/rys1_kaskada.svg
"""
import math
import os

W, H = 1040, 560
BLUE, ORANGE, RED, GREY, INK = "#1f5fa8", "#c8641a", "#b3261e", "#9aa3ad", "#1b1f24"

# (a) diagram Grotriana, skala energii schematyczna
YA = {1: 470, 2: 335, 3: 220, 4: 105}
XA = {0: 95, 1: 215, 2: 335, 3: 455}
HALF = 34
HALO = 'stroke="#ffffff" stroke-width="3.5" paint-order="stroke"'
# polozenie podpisu na strzalce: (t wzdluz strzalki, dx, dy, kotwica)
LABPOS = {((3, 1), (2, 0)): (0.55, -8, 0, "end"),
          ((4, 1), (2, 0)): (0.30, -8, 0, "end"), ((4, 2), (3, 1)): (0.45, 8, -4, "start")}
# (n, l) -> (n', l', udzial QM, udzial silnika lub None)
CHANNELS = [
    ((2, 1), (1, 0), 1.000, None),
    ((3, 0), (2, 1), 1.000, None),
    ((3, 1), (1, 0), 0.882, 0.882),
    ((3, 1), (2, 0), 0.118, 0.118),
    ((3, 2), (2, 1), 1.000, None),
    ((4, 0), (2, 1), 0.584, None),
    ((4, 0), (3, 1), 0.416, None),
    ((4, 1), (1, 0), 0.839, 0.841),
    ((4, 1), (2, 0), 0.119, 0.122),
    ((4, 1), (3, 0), 0.038, 0.038),
    ((4, 1), (3, 2), 0.004, None),
    ((4, 2), (2, 1), 0.746, 0.734),
    ((4, 2), (3, 1), 0.254, 0.266),
    ((4, 3), (3, 2), 1.000, None),
]
LABELLED = {((3, 1), (2, 0)), ((4, 1), (2, 0)), ((4, 2), (3, 1))}
RATIO = {(2, 1): "1,050", (3, 1): "1,060", (4, 1): "1,061", (3, 0): "0,977", (4, 0): "1,016"}

# (b) plaszczyzna dzialan: L/hbar poziomo, n = J/hbar pionowo
X0, SX = 640, 78
YB = {1: 470, 2: 350, 3: 230, 4: 110}


def xb(L):
    return X0 + SX * L


out = []
add = out.append
add('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" width="%d" height="%d" '
    'font-family="Georgia, \'Times New Roman\', serif" font-size="15" fill="%s">' % (W, H, W, H, INK))
add('<rect width="100%" height="100%" fill="#ffffff"/>')
add('<defs>')
for name, col in (("b", BLUE), ("o", ORANGE), ("r", RED)):
    add('<marker id="h%s" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="11" markerHeight="11" '
        'orient="auto-start-reverse" markerUnits="userSpaceOnUse"><path d="M0,0 L10,5 L0,10 z" fill="%s"/></marker>' % (name, col))
add('</defs>')

# ---------- panel (a)
late = []
add('<text x="20" y="34" font-size="18" font-weight="bold">(a)</text>')
add('<text x="52" y="34" font-size="16">kaskada E1 wodoru, n ≤ 4</text>')
for l, x in XA.items():
    add('<text x="%d" y="62" text-anchor="middle" font-style="italic" font-size="17">%s</text>' % (x, "spdf"[l]))
for n, y in YA.items():
    add('<text x="22" y="%d" font-size="15">n = %d</text>' % (y + 5, n))
    for l in range(n):
        add('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="%s" stroke-width="2.4"/>'
            % (XA[l] - HALF, y, XA[l] + HALF, y, INK))
        if (n, l) in RATIO:
            late.append('<text x="%d" y="%d" font-size="11.5" fill="#5b6570" %s>×%s</text>'
                % (XA[l] + HALF + 3, y + 4, HALO, RATIO[(n, l)]))


# przesuniecia koncow strzalek, zeby kanaly z jednego poziomu sie nie nakladaly
out_count = {}
for (i, f, *_r) in CHANNELS:
    out_count[i] = out_count.get(i, 0) + 1
used = {}
for (i, f, qm, eng) in CHANNELS:
    k = used.get(i, 0)
    used[i] = k + 1
    m = out_count[i]
    off = (k - (m - 1) / 2) * 14
    down = f[1] == i[1] - 1
    col, mk = (BLUE, "b") if down else (ORANGE, "o")
    x1, y1 = XA[i[1]] + off, YA[i[0]] + 3
    x2, y2 = XA[f[1]] + 0.5 * off, YA[f[0]] - 5
    wdt = 1.0 + 3.5 * qm
    add('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="%s" stroke-width="%.2f" '
        'stroke-opacity="0.9" marker-end="url(#h%s)"/>' % (x1, y1, x2, y2, col, wdt, mk))
    if (i, f) in LABELLED:
        t, ddx, ddy, anc = LABPOS[(i, f)]
        lx, ly = x1 + t * (x2 - x1), y1 + t * (y2 - y1)
        add('<text x="%.0f" y="%.0f" font-size="12.5" fill="%s" text-anchor="%s" %s>'
            '%s | %s</text>' % (lx + ddx, ly + ddy, col, anc, HALO, ("%.3f" % qm).replace(".", ","),
                                ("%.3f" % eng).replace(".", ",")))

out.extend(late)
add('<text x="20" y="510" font-size="12.5">grubość ∝ udział QM; liczby: udział QM | silnik; ×: czas życia silnik/QM</text>')
add('<line x1="22" y1="533" x2="56" y2="533" stroke="%s" stroke-width="3" marker-end="url(#hb)"/>' % BLUE)
add('<text x="64" y="538" font-size="12.5">Δl = −1, mod współbieżny</text>')
add('<line x1="250" y1="533" x2="284" y2="533" stroke="%s" stroke-width="3" marker-end="url(#ho)"/>' % ORANGE)
add('<text x="292" y="538" font-size="12.5">Δl = +1, mod przeciwbieżny</text>')

# ---------- panel (b)
add('<text x="560" y="34" font-size="18" font-weight="bold">(b)</text>')
add('<text x="592" y="34" font-size="16">płaszczyzna działań (L, J)</text>')
add('<line x1="%d" y1="500" x2="%d" y2="500" stroke="%s" stroke-width="1.2"/>' % (X0, xb(4.9), INK))
add('<line x1="%d" y1="500" x2="%d" y2="72" stroke="%s" stroke-width="1.2"/>' % (X0, X0, INK))
for L in range(0, 5):
    add('<line x1="%d" y1="500" x2="%d" y2="506" stroke="%s"/>' % (xb(L), xb(L), INK))
    add('<text x="%d" y="523" text-anchor="middle" font-size="13">%d</text>' % (xb(L), L))
add('<text x="%d" y="546" text-anchor="middle" font-size="14">|L| / ħ</text>' % xb(2.45))
add('<text x="604" y="80" font-size="14">J/ħ</text>')
for n, y in YB.items():
    add('<text x="%d" y="%d" text-anchor="end" font-size="13">%d</text>' % (X0 - 8, y + 4, n))
    # przedzial istnienia elipsy 0 <= L <= n hbar
    add('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="%s" stroke-width="7" stroke-opacity="0.35" '
        'stroke-linecap="round"/>' % (xb(0), y, xb(n), y, GREY))
    for l in range(n):
        add('<circle cx="%.1f" cy="%d" r="5.2" fill="%s"/>' % (xb(l + 0.5), y, INK))
    add('<circle cx="%.1f" cy="%d" r="5.2" fill="#ffffff" stroke="%s" stroke-width="1.6"/>'
        % (xb(n + 0.5), y, RED))
# granica L = J
add('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-dasharray="5 4" stroke-width="1.3"/>'
    % (xb(0.7), YB[1] + 0.3 * 120, xb(4.45), YB[4] - 0.45 * 120, GREY))
add('<text x="%.0f" y="%.0f" font-size="13" fill="%s" text-anchor="end">L = J (koło)</text>' % (xb(4.05), YB[4] - 32, GREY))


def arrow_b(a, b, col, mk, dash=False, bend=0.0):
    (L1, n1), (L2, n2) = a, b
    x1, y1, x2, y2 = xb(L1), YB[n1], xb(L2), YB[n2]
    dx, dy = x2 - x1, y2 - y1
    d = math.hypot(dx, dy)
    x1, y1 = x1 + 8 * dx / d, y1 + 8 * dy / d
    x2, y2 = x2 - 9 * dx / d, y2 - 9 * dy / d
    cx, cy = (x1 + x2) / 2 - bend * dy / d, (y1 + y2) / 2 + bend * dx / d
    add('<path d="M%.1f,%.1f Q%.1f,%.1f %.1f,%.1f" fill="none" stroke="%s" stroke-width="2.4"%s '
        'marker-end="url(#h%s)"/>' % (x1, y1, cx, cy, x2, y2, col,
                                     ' stroke-dasharray="6 4"' if dash else "", mk))
    t = 0.68
    return ((1 - t) ** 2 * x1 + 2 * t * (1 - t) * cx + t * t * x2,
            (1 - t) ** 2 * y1 + 2 * t * (1 - t) * cy + t * t * y2)


for a, b in (((3.5, 4), (2.5, 3)), ((2.5, 3), (1.5, 2)), ((1.5, 2), (0.5, 1))):
    arrow_b(a, b, BLUE, "b")
arrow_b((1.5, 4), (0.5, 2), BLUE, "b", bend=-14)
arrow_b((0.5, 3), (1.5, 2), ORANGE, "o", bend=10)
cx, cy = arrow_b((2.5, 3), (1.5, 1), RED, "r", dash=True, bend=-18)
add('<text x="%.0f" y="%.0f" font-size="22" font-weight="bold" fill="%s" text-anchor="middle">×</text>' % (cx, cy + 7, RED))
add('<text x="%.0f" y="%.0f" font-size="12.5" fill="%s">3d → 1p: L′ = 1,5ħ &gt; n′ħ,</text>'
    % (xb(1.95), YB[1] - 22, RED))
add('<text x="%.0f" y="%.0f" font-size="12.5" fill="%s">elipsa nie istnieje</text>'
    % (xb(1.95), YB[1] - 7, RED))
add('<text x="%.0f" y="%.0f" font-size="12.5" fill="%s">yrast 4f → 3d → 2p → 1s</text>'
    % (xb(2.75), YB[3] + 52, BLUE))
add('<text x="%.0f" y="%.0f" font-size="12.5" fill="%s" %s>4p → 2s (Δn = 2)</text>' % (xb(1.5) + 6, YB[4] + 72, BLUE, HALO))
add('<text x="%.0f" y="%.0f" font-size="12.5" fill="%s">3s → 2p</text>' % (xb(0.05), YB[2] + 58, ORANGE))
add('<circle cx="%d" cy="531" r="5" fill="%s"/>' % (xb(2.8), INK))
add('<text x="%d" y="536" font-size="12.5">L = (l+½)ħ</text>' % (xb(2.8) + 10))
add('<circle cx="%d" cy="531" r="5" fill="#ffffff" stroke="%s" stroke-width="1.6"/>' % (xb(4.1), RED))
add('<text x="%d" y="536" font-size="12.5">l = n</text>' % (xb(4.1) + 10))
add('</svg>')

path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "rys1_kaskada.svg")
with open(path, "w") as fh:
    fh.write("\n".join(out) + "\n")
print(path)
