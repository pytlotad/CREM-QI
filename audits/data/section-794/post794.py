# Audit 794 POST HOC: (a) Omega change from the start fix vs the energy-offset prediction -1,5 x (E0/E_S - 1) (hydrogen seed 5, A -> B);
# (b) residual deficit after the fix vs energy drift over 200 orbits (B trace).  usage: python3 post794.py
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()): import ana794 as A
pr = A.pr
pr('=== (a) Omega(B)/Omega(A) - 1 wobec -1,5 x (E0/E_S - 1), wodor ziarno 5 ===')
for n, L in A.TORI:
    p = lambda c: f'{A.O}/{c}/proton_n{n}L{A.lab(L)}_5.err'
    rs = A.relstart(p('B')); d = A.omega(p('B'))[0]/A.omega(p('A'))[0] - 1; pred = -1.5*(rs['E0']/rs['ES'] - 1)
    pr(f'  n {n} L {L:3} (e {A.ecc(n, L):.3f}): zmiana {d:+.4%} | przewidywanie {pred:+.4%} | roznica {d - pred:+.4%}')
pr('=== (b) resztkowy deficyt B (srednie 6 biegow, wobec n3L3) i dryf energii B (wodor ziarno 5) ===')
KM = A.RES['B']['KM'] if hasattr(A, 'RES') else None
for n, L in A.TORI:
    p = f'{A.O}/B/proton_n{n}L{A.lab(L)}_5.err'; tb = A.trace(p); rs = A.relstart(p)
    pr(f'  n {n} L {L:3} (e {A.ecc(n, L):.3f}): deficyt {KM[(n, L)]/KM[(3, 3)] - 1:+.4%} | dryf E(200)/E(1) - 1 {tb[-1]["E"]/tb[0]["E"] - 1:+.2e} | '
       f'1,5 x polowa dryfu {0.75*(tb[-1]["E"]/tb[0]["E"] - 1):+.4%}')
